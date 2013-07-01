#include "codonUtils.h"
#include "numRec.h"
#include <algorithm>



void printHelp(){

		cout <<"+-------------------------------------------------------+"<<endl;
		cout <<"Input:"<<endl;
		cout <<"---------------------------------------------------------"<<endl;
		cout <<"-i    input codon-aligned sequence file "<<endl;
		cout <<"      (accepted formats:  Fasta, Clustal)"<<endl;
		cout <<"-q    name of query sequence (default=1st in file)"<<endl;
		cout <<"---------------------------------------------------------"<<endl;
		cout <<"Advanced options:"<<endl;
		cout <<"---------------------------------------------------------"<<endl;
		cout <<"-u    input user tree in Newick format"<<endl;				
		cout <<"-g    genetic code (default: nuc. standard)"<<endl;
		cout <<"	Nuclear:"<<endl;
		cout <<"	0: standard , 1:  Blepharisma"<<endl;
		cout <<"	2:  Ciliate, 3:  Euplotid"<<endl;
		cout <<"	Mitochondria:"<<endl;
		cout <<"	4: Vertebrate, 5: Ascidian, 6: Echinoderm"<<endl;
		cout <<"	7: Flatworm, 8: Invertebrate"<<endl;
		cout <<"	9: Protozoan, 10: Yeast"<<endl;
		cout <<"-m    method:"<<endl;
		cout <<"	-mb    bayesian	(default)"<<endl;
		cout <<"	-ml    maximum likelihood"<<endl;
		cout <<"-d    prior bayesian distribution"<<endl;
		cout <<"	-db    beta+w (default)"<<endl;
		cout <<"	-dg    gamma"<<endl;
		cout <<"-n  No. of categories for discrete distr.(default=8)"<<endl;
		cout <<"-bn   no branch length optimization"<<endl;
		cout <<"-e    epsilon for likelihood optimization (default=0.1)"<<endl;
		cout <<"      (the smaller the value, the higher the precision)"<<endl;
		cout <<"-j    number of optimization iterations (default=5)"<<endl;
		cout <<"--------------------------------------------------------"<<endl;
		cout <<"For fixing the parameters, or running a specific model:"<<endl;
		cout <<"--------------------------------------------------------"<<endl;
		cout <<"(Note that if one of the below options is not used,"<<endl;
		cout <<"the default run will be of the M8 model)"<<endl;
		cout <<"-w    initial value of additional w category (-w1=M8a model)"<<endl;
		cout <<"-p    initial probability of beta distribution"<<endl;
		cout <<"-a    initial alpha value"<<endl;
		cout <<"-x    initial beta value"<<endl;
		cout <<"-k    initial kappa value"<<endl;
		cout <<"-Fw   fixed value of additional omega value"<<endl;
		cout <<"-Fp   fixed probability of beta distribution"<<endl;
		cout <<"-Fa   fixed alpha value"<<endl;
		cout <<"-Fx   fixed beta value"<<endl;
		cout <<"-Fk   fixed kappa"<<endl;
		cout <<"** For the M8a model, type -w1 -Fw"<<endl;
		cout <<"** For the M7 model, type -p1 -Fp"<<endl;
		cout <<"--------------------------------------------------------"<<endl;
		cout <<"Output files:"<<endl;
		cout <<"--------------------------------------------------------"<<endl;
		cout <<"-l    log file"<<endl;
		cout <<"-r    results output file with Ka/Ks, CI and posterior"<<endl;
		cout <<"-o    output file with likelihood and optimized params"<<endl;
		cout <<"-s    Rasmol script for coloring a 3D molecule, if available"<<endl;
		cout <<"-c    output color bin file (site to color, according to web server colors)"<<endl;
		cout <<"-t    output tree file"<<endl;
		cout <<"---------------------------------------------------------"<<endl;
		cout <<"-h or -? or -H     help"<<endl;
		cout <<"lowercase and uppercase letters are both ok"<<endl;
		cout <<"---------------------------------------------------------+"<<endl;
}


//check that the input sequences are divisable by 3
void checkInputSeqLength(string codonFile){
	nucleotide alph;
	ifstream in(codonFile.c_str());
	sequenceContainer inputSc = recognizeFormat::readUnAligned(in, &alph);
	in.close();
	int i;
    for (i = 0; i < inputSc.numberOfSeqs(); ++i){
		int seqLen = inputSc[i].seqLen();
		if ((seqLen % 3) != 0){
			string textToPrint = "USER ERROR: unable to read sequence: " + inputSc[i].name() + "\nSequence length is not divisable by three";
			errorMsg::reportError(textToPrint);
		}
	}
}

//this function convert codon sequences to amino sequences.
sequenceContainer convertCodonToAmino(sequenceContainer &codonSc,codon *codonAlph){
	amino aaAlph;
	sequenceContainer aaSc;
	for (int i = 0; i < codonSc.numberOfSeqs(); ++i){
		sequence codonSeq = codonSc[i];
		sequence aaSeq("", codonSeq.name(), codonSeq .remark(), codonSeq.id(), &aaAlph);
		for (int pos = 0; pos < codonSeq .seqLen(); ++pos)
            aaSeq.push_back(codonUtility::aaOf(codonSeq[pos],*codonAlph));
		aaSc.add(aaSeq);
	}
	if (codonSc.numberOfSeqs() != aaSc.numberOfSeqs())
		errorMsg::reportError("RevTrans: number of codon and Amino sequences is not the same");
	
	return aaSc;
}

// normalize the Q matrix so average rate of substitution = 1
void normalizeMatrices(vector<stochasticProcess> & spVec,const distribution * forceDistr){
	MDOUBLE sumPijQij=0.0;
	int categor;
	for ( categor=0; categor<forceDistr->categories();categor++)
		sumPijQij+=forceDistr->ratesProb(categor)*static_cast<wYangModel*>(spVec[categor].getPijAccelerator()->getReplacementModel())->sumPijQij();	
	if (sumPijQij ==0){
		errorMsg::reportError("Error in normalizeMatrices - sumPijQij=0");
	}
	for (categor=0; categor<forceDistr->categories();categor++)
		static_cast<wYangModel*>(spVec[categor].getPijAccelerator()->getReplacementModel())->norm(1/sumPijQij);

}

Vdouble freqCodonF3x4(const sequenceContainer &nucSc, codon * coAlph){
	VVdouble nucFeqPos(3);
	int pos= 0;
	int nPos = 0;
	for (nPos=0;nPos<3;nPos++)
		nucFeqPos[nPos].resize(nucSc.alphabetSize(),0.0);

	sequenceContainer::constTaxaIterator tIt;
	sequenceContainer::constTaxaIterator tItEnd;
	tIt.begin(nucSc);
	tItEnd.end(nucSc);
	while (tIt!= tItEnd) {
		pos = 0;
		sequence::constIterator sIt;
		sequence::constIterator sItEnd;
		sIt.begin(*tIt);
		sItEnd.end(*tIt);
		while (sIt != sItEnd) {
			if ((*sIt >= 0) && (*sIt <nucFeqPos[pos%3].size())) ++nucFeqPos[pos%3][(*sIt)];
			if (*sIt == 4) ++nucFeqPos[pos%3][3]; //for T (4) to U (3)
			++sIt;
			++pos;
		}
		++tIt;
	}
	for (nPos=0;nPos<3;nPos++)
		changeCountsToFreqs(nucFeqPos[nPos]);


	Vdouble freqCodon(coAlph->size(),0.0);

	nucleotide n;
	for (int c = 0; c<freqCodon.size();c++){

		string s = coAlph->fromInt(c);
		int nuc0 = n.fromChar(s[0]);
		int nuc1 = n.fromChar(s[1]);
		int nuc2 = n.fromChar(s[2]);
		freqCodon[c] = nucFeqPos[0][nuc0]*nucFeqPos[1][nuc1]*nucFeqPos[2][nuc2];
	}

	MDOUBLE sum=0;
	for (int i=0;i<coAlph->size();i++){
		sum+=freqCodon[i];
	}
	MDOUBLE stopFreq  = 1.0 - sum;
	MDOUBLE  ep = stopFreq/coAlph->size();
	for (int i=0;i<coAlph->size();i++){
		freqCodon[i]+=ep;
	}

	return freqCodon;


}


/***********************************************
  The following functions are useful for the selecton server, for creating a 
  Rasmol script and for setting the color value of each site
 ***********************************************/


// Positive significant in color dark yellow, non-sig. positive selection - light yellow. 
// Purifying selection in shades of bordeaux
vector<vector<int> > create7ColorValues(){
	vector<vector<int> > colorsValue;
	colorsValue.resize(7);
	for (int i=0;i<7;i++)
		colorsValue[i].resize(3);
	// RGB values of the differnt color bins
	colorsValue[0][0] = 255; //yellow positive significant
	colorsValue[0][1] = 220 ;
	colorsValue[0][2] = 0;

	colorsValue[1][0] =255 ; //light yellow - not  significant positive selection
	colorsValue[1][1] = 255;
	colorsValue[1][2] = 120;
	
	//three categories of not significant negative selection according to bordeaux shades (colors like conseq/consurf)

	colorsValue[2][0] = 255; //white  
	colorsValue[2][1] = 255;
	colorsValue[2][2] = 255;

	colorsValue[3][0] = 252;
	colorsValue[3][1] = 237;
	colorsValue[3][2] = 244;

	colorsValue[4][0] = 250;
	colorsValue[4][1] = 201;
	colorsValue[4][2] = 222;

	colorsValue[5][0] = 240;
	colorsValue[5][1] = 125;
	colorsValue[5][2] = 171;
	
	//significant negative selection
	colorsValue[6][0] = 130; 
	colorsValue[6][1] = 67;
	colorsValue[6][2] = 96;

	return colorsValue;
}

//this functions creates a rasmol script (assumes positions are the same between the alignment and the PDB)
void outToRasmolFile(string fileName,vector<int>& color4Site){
	ofstream out(fileName.c_str());	
	vector<vector<int> > colorsValue = create7ColorValues();
	int numberOfColor = colorsValue.size();
	vector<vector<int> > colors; //for each color (1-9/3) holds vector of sites.
	colors.resize(numberOfColor+1);
	int i;
	for (i=0;i<color4Site.size();i++){
		int color=color4Site[i];
		if (color>numberOfColor){
			errorMsg::reportError("Error in outToColorFile - unknown color");
		}
		colors[color].push_back(i+1); //add site (position in the vector +1)	
	}
	out<<"select all"<<endl;
	out<<"color [200,200,200]"<<endl<<endl;
	
	for (int c=1;c<numberOfColor+1;c++){
		out<<"select ";
		for (i=0;i<colors[c].size();i++){
			if (i==0)
				out<<colors[c][i];
			else if ((i+1)%6==0)
				out<<endl<<"select selected or "<<colors[c][i];
			 
			else out<<" , "<<colors[c][i];
		}
		out<<endl<<"select selected and :a"<<endl;
		out<<"color [" <<colorsValue[c-1][0]<<","<<colorsValue[c-1][1]<<","<<colorsValue[c-1][2]<<"]"<<endl;
		out<<"spacefill"<<endl<<endl;
	}
	
	out.close();
}


// a file with color-coding from Ka/Ks values to color-bins
void kaks2Color(const Vdouble & kaksVec, const Vdouble &lowerBoundV,
				const sequence & refSeq, string fileName,codon *co) {
	vector<int> colors;
	int numOfSitesinAln = kaksVec.size();
	Vdouble negativesKaksVec,negativesSite;
	negativesKaksVec.clear();
	negativesSite.clear();
	int i,gapsInRefSeq=0;

	for (i=0;i<numOfSitesinAln;i++){
		if (codonUtility::aaOf(refSeq[i],*co) == -1) gapsInRefSeq++; 
	}

	// first dealing with positive selection
	colors.resize(numOfSitesinAln-gapsInRefSeq);
	int gap=0;
	for (i=0;i<numOfSitesinAln;i++){
		if (codonUtility::aaOf(refSeq[i],*co) == -1){
			gap++;
			continue;
		}
		if (lowerBoundV[i]>1) // color 1 (positive selection) : if confidence interval lower bound > 1
			colors[i-gap]=1;
		else if (kaksVec[i]>1) // color 2(positive selection) : "non-significant"
			colors[i-gap]=2;
		else  {
			negativesKaksVec.push_back(kaksVec[i]);  //add the value of kaks < 1
			negativesSite.push_back(i-gap);   //add the number of site of the kaks 
		}
	
	}

	// now dealing with purifying selection
	Vdouble orderVec = negativesKaksVec;
	if (orderVec.size()>0) // this is since once the whole protein was positive selection... (anomaly)
		sort(orderVec.begin(), orderVec.end());  //sort the kaks values to be divided to 5 groups
	MDOUBLE percentileNum = 5.0;
	int percentileNumInt = 5;
	Vdouble maxScoreForPercentile(percentileNumInt);
	if (orderVec.size()>0) {
		maxScoreForPercentile[0] = orderVec[0]; 
		for (int c = 1; c < percentileNumInt; ++c){
			int place = (int)((c / percentileNum) * negativesKaksVec.size());
			MDOUBLE maxScore = orderVec[place];
			maxScoreForPercentile[c] = maxScore;
		}
	}

	//loop over all the Ka/Ks < 1  
	for (int j=0; j < negativesKaksVec.size(); ++j){
			MDOUBLE r = negativesKaksVec[j]; //the kaks of the site.
			int s = (int)negativesSite[j];  //the  site.
			if (r > maxScoreForPercentile[4]) 
					colors[s] = 3;
			else if (r > maxScoreForPercentile[3]) 
					colors[s] = 4;
			else if (r> maxScoreForPercentile[2])
					colors[s] = 5;
			else if (r > maxScoreForPercentile[1])
					colors[s] = 6;
			else if (r >= maxScoreForPercentile[0])
					colors[s] = 7;
	}
	//print to file
	ofstream out(fileName.c_str());
	gap=0;
	amino aminoAcid;
	LOG(5,<<"Printing selection color bins to file"<<endl);
	for (i=0;i<refSeq.seqLen();i++){	 
		int aa = codonUtility::aaOf(refSeq[i], *co);
		if (aa==-1){
			gap++;
			continue;
		}
		string aaStr = aminoAcid.fromInt(aa);
		out<<i+1-gap <<"\t"<<aaStr<<"\t"<<colors[i-gap];
		out<<endl;
	}
	out.close();
}
