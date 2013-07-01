#include "mainbb.h"

#include "aaJC.h"
#include "amino.h"
#include "bbAlg.h"
#include "bestAlpha.h"
#include "bblEM.h"
#include "chebyshevAccelerator.h"
#include "clustalFormat.h"
#include "computeMarginalReconstruction.h"
#include "distanceTable.h"
#include "fastaFormat.h"
#include "gammaDistribution.h"
#include "jointNoGamma.h"
#include "likeDist.h"
#include "logFile.h"
#include "maseFormat.h"
#include "molphyFormat.h"
#include "nexusFormat.h"
#include "nucleotide.h"
#include "nucJC.h"
#include "gtrModel.h"
#include "nj.h"
#include "phylipFormat.h"
#include "readDatMatrix.h"
#include "recognizeFormat.h"
#include "trivialAccelerator.h"
#include "uniDistribution.h"

//For the codon part
#include "bestAlphaAndK.h"
#include "codonUtils.h"


#include <fstream>
#include <iostream>
using namespace std;

mainbb::mainbb(int argc, char* argv[]) {
	fillOptionsParameters(argc,argv);
	myLog::setLog(_options->reportFile,10);
	printBBProjectInfo();
	printSearchParameters();
	getStartingSequenceData();
	getStartingStochasticProcess();
	getStartingEvolTreeTopology();
	_et.rootToUnrootedTree();
	//_et.createFlatLengthMatrix(0.001); // TO BE USED FOR TESTING ONLY.
	if (_options->modelName == bb_options::nyCodon)		
		getStartingBLAndModelParam(); //for NY codon Models 
	else 
		getStartingBranchLengthsAndAlpha();
	printOutputTree();
	if (_options->doJoint) {
		if (_options->distributionName == bb_options::gam) {
            findAncestralSequencesGammaJoint();
        } else {
			findAncestralSequencesHomJoint();		
		}
	}
	else{
		getMarginalReconstruction();
	}
	myLog::endLog();
}

void mainbb::printAncestralSequencesGammaJoint() {
	replaceSequences(_resulutingJointReconstruction,_originSc);
	ofstream out(_options->outFile_seq_joint.c_str());
	//out<<"sequences of the joint reconstruction, model: "<<_options->modelNameStr()<<endl;
	switch (_options->seqOutputFormat){
		case (bb_options::mase)   : maseFormat::write(out,_resulutingJointReconstruction); break;
		case (bb_options::fasta)  : fastaFormat::write(out,_resulutingJointReconstruction); break;
		case (bb_options::clustal): clustalFormat::write(out,_resulutingJointReconstruction); break;
		case (bb_options::phylip) : phylipFormat::write(out,_resulutingJointReconstruction); break;
		case (bb_options::molphy) : molphyFormat::write(out,_resulutingJointReconstruction); break;
		case (bb_options::nexus)  : nexusFormat::write(out,_resulutingJointReconstruction); break;
	}
	out.close();
}

mainbb::~mainbb() {
	if (_alph) delete _alph;
	if (_options) delete _options;
}

void mainbb::getStartingEvolTreeTopology(){
	if (_options->treefile=="") {
		getStartingNJtreeNjMLdis();
	}
	else getStartingTreeFromTreeFile();
}



void mainbb::getStartingNJtreeNjMLdis() {
	// note that here ALWAYS, the ML distances are computed using
	// an homogenous rate distribution.
	uniDistribution lUni;
//	const pijAccelerator* lpijAcc = _sp->getPijAccelerator();// note this is just a copy of the pointer.
	const pijAccelerator* lpijAcc = _spVec[0].getPijAccelerator();// note this is just a copy of the pointer.
	stochasticProcess lsp(&lUni,lpijAcc);

	likeDist pd1(lsp,0.01);
	VVdouble disTab;
	vector<string> vNames;
	giveDistanceTable(&pd1,
					   _sc,
					   disTab,
					   vNames);
	getStartingTreeNJ_fromDistances(disTab,vNames);
} 

void mainbb::getStartingTreeNJ_fromDistances(const VVdouble& disTab,
	const vector<string>& vNames) {
	NJalg nj1;
	_et= nj1.computeTree(disTab,vNames);

}
	
void mainbb::getStartingTreeFromTreeFile(){
	_et= tree(_options->treefile);
	if (!_et.withBranchLength()) {
		_et.createFlatLengthMatrix(0.05);
		_options->optimizeBrLenOnStartingTree = true;
	}
}

void mainbb::getStartingBranchLengthsAndAlpha(){
	if (_options->distributionName == bb_options::hom) {
		if (_options->optimizeBrLenOnStartingTree == true) {
			cout<<"Optimizing branch lengths (Homogenuos model)..."<<endl;
			bblEM bblem1(_et,_sc,_spVec[0],NULL);
			//bblEM bblem1(_et,_sc,*_sp,NULL);
			//brLenOptEM::optimizeBranchLength1G_EM(_et,_sc,*_sp,NULL);
		}
	}
	else { // GAMMA MODEL!
		// Here we want to optimize branch lengths with a gamma model.
		// there are three options:
		//(1) User provides the alpha and no bbl.
		//(2) User provides the alpha and bbl
		//(3) Alpha is optimized from the data and bbl.


		// option 1 will not enter to any of these options.
		if ((_options->userProvideAlpha == true) && (_options->optimizeBrLenOnStartingTree == true)) {
			cout<<"Optimizing branch lengths (Gamma model, user alpha)..."<<endl;
			MDOUBLE intitalAlpha = 1.0;
			static_cast<gammaDistribution*>(_spVec[0].distr())->setAlpha(intitalAlpha);
			bblEM bblem1(_et,_sc,_spVec[0],NULL);
			//static_cast<gammaDistribution*>(_sp->distr())->setAlpha(intitalAlpha);
			//bblEM bblem1(_et,_sc,*_sp,NULL);
			//brLenOptEM::optimizeBranchLength1G_EM(_et,_sc,*_sp,NULL);
		} 
		else if ((_options->userProvideAlpha == true) && (_options->optimizeBrLenOnStartingTree == false)) {
			return;
		}
		else if (_options->userProvideAlpha == false)  {
			cout<<"Optimizing branch lengths and alpha (Gamma model) ..."<<endl;
			bestAlphaAndBBL bbl2(_et,_sc,_spVec[0]);
		}
	}
}

void mainbb::getStartingStochasticProcess() {
	int numberOfCategories = _options->gammaCategies;
	MDOUBLE alpha = _options->gammaPar;
	if (_options->distributionName == bb_options::hom) {
		numberOfCategories = 1; // forcing homogenous model.
		alpha = 1.0;
		cout<<"Using homogenous model (no among site rate variation)"<<endl;
	} else {
		cout<<"Using a Gamma model with: "<<numberOfCategories<<" discrete categories "<<endl;
	}
	distribution *dist = new gammaDistribution(alpha,numberOfCategories);
	replacementModel *probMod=NULL;
	pijAccelerator *pijAcc=NULL;
	switch (_options->modelName){
		case (bb_options::day): 
			probMod=new pupAll(datMatrixHolder::dayhoff);
			if (_options->useChebyshev == true) {
				pijAcc = new chebyshevAccelerator(probMod);
			} else {
				pijAcc = new trivialAccelerator(probMod);
			}
			cout<<"Amino acid replacement matrix is Dayhoff"<<endl;
			break;
		case (bb_options::jtt):
			probMod=new pupAll(datMatrixHolder::jones);
			if (_options->useChebyshev == true) {
				pijAcc = new chebyshevAccelerator(probMod);
			} else {
				pijAcc = new trivialAccelerator(probMod);
			}
			cout<<"Amino acid replacement matrix is JTT"<<endl;
			break;
		case (bb_options::lg):
			probMod=new pupAll(datMatrixHolder::lg);
			if (_options->useChebyshev == true) {
				pijAcc = new chebyshevAccelerator(probMod);
			} else {
				pijAcc = new trivialAccelerator(probMod);
			}
			cout<<"Amino acid replacement matrix is LG"<<endl;
			break;
		case (bb_options::rev):
			probMod=new pupAll(datMatrixHolder::mtREV24);
			if (_options->useChebyshev == true) {
				pijAcc = new chebyshevAccelerator(probMod);
			} else {
				pijAcc = new trivialAccelerator(probMod);
			}
			cout<<"Amino acid replacement matrix is mtREV24"<<endl;
			 break;
		case (bb_options::wag):
			probMod=new pupAll(datMatrixHolder::wag);
			if (_options->useChebyshev == true) {
				pijAcc = new chebyshevAccelerator(probMod);
			} else {
				pijAcc = new trivialAccelerator(probMod);
			}
			cout<<"Amino acid replacement matrix is WAG"<<endl;
			break;
		case (bb_options::cprev):
			probMod=new pupAll(datMatrixHolder::cpREV45);
			if (_options->useChebyshev == true) {
				pijAcc = new chebyshevAccelerator(probMod);
			} else {
				pijAcc = new trivialAccelerator(probMod);
			}
			cout<<"Amino acid replacement matrix is cpREV45"<<endl;
			break;
		case (bb_options::empiriCodon):
			probMod=new pupAll(datMatrixHolder::empiriCodon,61);
			if (_options->useChebyshev == true) {
				pijAcc = new chebyshevAccelerator(probMod,61);
			} else {
				pijAcc = new trivialAccelerator(probMod);
			}
			cout<<"Codon replacement matrix is empiriCodon of adrian"<<endl;
			break;
		case (bb_options::nucjc):
			probMod=new nucJC;
			pijAcc = new trivialAccelerator(probMod);
			cout<<"Nucleotide substitution model is Jukes and Cantor"<<endl;
			break;
		case (bb_options::nucgtr):
			{
			nucleotide nucAlph;
			Vdouble freq = computeGTRFreq(nucAlph);
			probMod=new gtrModel(freq, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25);
			pijAcc = new trivialAccelerator(probMod);
			cout<<"Nucleotide substitution model is General time Reversible"<<endl;
			}
			break;
		case (bb_options::aajc):
			probMod=new aaJC; pijAcc = new trivialAccelerator(probMod);
			cout<<"Amino acid replacement matrix is Jukes and Cantor"<<endl;
			break;
				//this part for the codon model c & w init as with no selection	
		case (bb_options::nyCodon):
			{
				codon codonAlph;
				Vdouble freq = computeFreq(codonAlph);
				probMod = new wYangModel(1.0,1.0,freq, 0, &codonAlph); 
				pijAcc = new trivialAccelerator(probMod);
				cout<<"Codon replacement matrix is NY model"<<endl;
			}
			break;
		default:
			errorMsg::reportError("this probablistic model is not yet available");
	}
	stochasticProcess sp(dist, pijAcc);
	_spVec.push_back(sp);
	if (probMod) delete probMod;
	if (pijAcc) delete pijAcc;
	if (dist) delete dist;
}

void mainbb::printOutputTree() {
	ofstream f;
	string fileName1=_options->outTreeFileNewick;
	f.open(fileName1.c_str());
	_et.output(f,tree::PHYLIP,true);
	//_et.output(f,tree::PHYLIP,false);
	f.close();
	cout<<"The tree in 'Newick tree format' (with the internal nodes labeled)\nwas written to a file name called "<<fileName1<<endl;
	fileName1 = _options->outTreeFileAncestor;
	f.open(fileName1.c_str());
	_et.output(f,tree::ANCESTOR);
	f.close();
	cout<<"The tree in 'ANCESTOR tree format' was written to a file name called "<<fileName1<<endl;
}

void mainbb::fillOptionsParameters(int argc, char* argv[]) {
	_options = new bb_options(argc, argv);
}

void mainbb::getStartingSequenceData(){	
 	if (_options->alphabet_size==4) _alph = new nucleotide;
	else if (_options->alphabet_size == 20) _alph = new amino;
	else if (_options->alphabet_size == 61) _alph = new codon;
	else errorMsg::reportError("no such alphabet in function rate4site::getStartingSequenceData");
   
	ifstream fstream1(_options->seqfile.c_str());
	_sc = recognizeFormat::read(fstream1,_alph);
	_originSc = _sc;
	_sc.changeGaps2MissingData();
}
	
void mainbb::printSearchParameters() {
	if (_options->verbose) {
		LOG(1,<<"\nBB parameters: "<<endl);
		LOG(1,<<endl);
		LOG(1,<<"-------------------------------------------------------------------------------"<<endl);
		LOG(1,<<endl);
		if (_options->treefile.size()>0) {LOG(1,<<"Tree file is: "<<_options->treefile<<endl)}
		else LOG(1,<<"Starting tree is the NJ tree "<<endl);
		if (_options->seqfile.size()>0) LOG(1,<<"Sequence file is: "<<_options->seqfile<<endl);
	}
}
             
void mainbb::printBBProjectInfo() {
	LOG(1,<<"*******************************************************************************"<<endl);
	LOG(1,<<"B&B: A Branch and Bound algorithm for Ancestral Sequence Reconstruction.       "<<endl);
	LOG(1,<<"For information, please send email to Tal Pupko: talp@post.tau.ac.il           "<<endl);
	LOG(1,<<"Ref: Pupko, T., Pe'er, I., Graur, D. Hasegawa, M., and Friedman N. 2002.       "<<endl);
	LOG(1,<<"A branch-and-bound algorithm for the inference of ancestral amino-acid         "<<endl);
	LOG(1,<<"sequences when the replacement rate varies among sites: Application to the     "<<endl);
	LOG(1,<<"evolution of five gene families. Bioinformatics 18: 1116-1123.                 "<<endl);
	LOG(1,<<"*******************************************************************************"<<endl);
	LOG(1,<<endl);
}

void mainbb::findAncestralSequencesGammaJoint() {
	bbAlg::boundMethod bm;
	if (_options->boundMethod == bb_options::max) bm=bbAlg::max;
	else if (_options->boundMethod == bb_options::sum) bm=bbAlg::sum;
	else if (_options->boundMethod == bb_options::both) bm=bbAlg::both;

	bbAlg bbAlg1(_et,_spVec,_sc,bm,_options->reportFile,_options->computeAgainExactTreshold,_forceDistr);
	cout<<"after bbAlg in findAncestralSequencesGammaJoint()"<<endl;
	//bbAlg bbAlg1(_et,*_sp,_sc,bm,_options->reportFile,_options->computeAgainExactTreshold);
	MDOUBLE res = bbAlg1.bbReconstructAllPositions(_resulutingJointReconstruction);
	cout<<" the likelihood of this reconstruction is: "<<res<<endl;
	bbAlg1.outputTheJointProbAtEachSite(_options->outFile_prob_joint);
	printAncestralSequencesGammaJoint();
}

void mainbb::findAncestralSequencesHomJoint() {
	//jointNoGamma jng(_et,*_sp,_sc);
	jointNoGamma jng(_et,_spVec[0],_sc);
	jng.compute();
	jng.outputTheJointProbAtEachSite(_options->outFile_prob_joint);
	sequenceContainer withAncestral = jng.getTheJointReconstruction();
	replaceSequences(withAncestral,_originSc);
	ofstream jointNoGammaReconstructionOutputFile(_options->outFile_seq_joint.c_str());
	//jointNoGammaReconstructionOutputFile<<"sequences of the joint reconstruction, model (hom): "<<_options->modelNameStr()<<endl;
	switch (_options->seqOutputFormat) {
	case bb_options::mase: 
			 maseFormat::write(jointNoGammaReconstructionOutputFile,withAncestral);
			break;
	case bb_options::molphy: 
			molphyFormat::write(jointNoGammaReconstructionOutputFile,withAncestral);
			break;
	case bb_options::clustal: 
			clustalFormat::write(jointNoGammaReconstructionOutputFile,withAncestral);
			break;
	case bb_options::fasta: 
			fastaFormat::write(jointNoGammaReconstructionOutputFile,withAncestral);
			break;
	case bb_options::phylip: 
			phylipFormat::write(jointNoGammaReconstructionOutputFile,withAncestral);
			break;
	case bb_options::nexus: 
			nexusFormat::write(jointNoGammaReconstructionOutputFile,withAncestral);
			break;
	default: errorMsg::reportError(" format not implemented yet in this version... ",1);
	}
}


void mainbb::getMarginalReconstruction(){
	//computeMarginalReconstruction cmr(_et,*_sp,_sc);
	computeMarginalReconstruction cmr(_et,_spVec,_sc);
	cmr.compute(_forceDistr);
	//cmr.compute();
	cmr.outputTheMarginalProbForEachCharForEachNode(_options->outFile_prob_marginal);
	sequenceContainer withAncestral = cmr.getResultingMarginalReconstruction();
	replaceSequences(withAncestral,_originSc);
	ofstream marginalReconstructionOutputFile(_options->outFile_seq_marginal.c_str());
	marginalReconstructionOutputFile<<"sequences of the marginal reconstruction, model: "<<_options->modelNameStr()<<endl;
	switch (_options->seqOutputFormat) {
	case bb_options::mase: 
			 maseFormat::write(marginalReconstructionOutputFile,withAncestral);
			break;
	case bb_options::molphy: 
			molphyFormat::write(marginalReconstructionOutputFile,withAncestral);
			break;
	case bb_options::clustal: 
			clustalFormat::write(marginalReconstructionOutputFile,withAncestral);
			break;
	case bb_options::fasta: 
			fastaFormat::write(marginalReconstructionOutputFile,withAncestral);
			break;
	case bb_options::phylip: 
			phylipFormat::write(marginalReconstructionOutputFile,withAncestral);
			break;
	case bb_options::nexus: 
			nexusFormat::write(marginalReconstructionOutputFile,withAncestral);
			break;
	default: errorMsg::reportError(" format not implemented yet in this version... ",1);
	}
	marginalReconstructionOutputFile.close();
}


//This part for NY codon model 
//for optomize the w yang model under gamma model and BBL
 void mainbb::getStartingBLAndModelParam()
 {
	 // GAMMA MODEL FOR W Yang Model
		// Here we want to optimize branch lengths with a gamma model.
		// there are three options:
		//(1) User provides the alpha and no bbl.
		//(2) User provides the alpha and bbl
		//(3) Alpha is optimized from the data and bbl.
		cout<<"Optimization of NY model with gamma - M5 in PAML"<<endl<<endl;
		createStochasticProcessVec();
		if ((_options->userProvideAlpha == true) && (_options->optimizeBrLenOnStartingTree == true)) {
			cout<<"Optimizing branch lengths & parametrs model: beta + k (Gamma model, user alpha)..."<<endl;			
			optimizeSelectonParameters bestParams(_et,_sc,_spVec,_forceDistr,true,true,false,false,false,true,false,3,3,0.01,0.01,0.1,20,20);
		}
		
		else if ((_options->userProvideAlpha == true) && (_options->optimizeBrLenOnStartingTree == false)) {
			cout<<"Optimizing parametrs model: k + beta (Gamma model, user alpha, user branch lengths)..."<<endl;
			optimizeSelectonParameters bestParams(_et,_sc,_spVec,_forceDistr,0,1,0,0,0,1,0);
		
		}
		else if (_options->userProvideAlpha == false)  {
			cout<<"Optimizing branch lengths and model parametrs alpha + beta +k (Gamma model) ... "<<endl;
			optimizeSelectonParameters bestParams(_et,_sc,_spVec,_forceDistr,1,1,0,0,0,0,0);
		}
 }


 void mainbb::createStochasticProcessVec()
 {
	 wYangModel * baseModel  = static_cast<wYangModel*>(_spVec[0].getPijAccelerator()->getReplacementModel());
	 wYangModel tmp(*baseModel);
	 _forceDistr = new generalGammaDistribution(_options->gammaPar,_options->gammaPar,_options->gammaCategies);
	_spVec.resize(_forceDistr->categories());
	uniDistribution dist;
	for (int categor=0; categor<_forceDistr->categories();categor++){
		wYangModel tmpModel(tmp);
		tmpModel.setW(_forceDistr->rates(categor));
		trivialAccelerator pijAcc(&tmpModel);
		stochasticProcess tmpSp(&dist,&pijAcc);
		_spVec[categor] = tmpSp;
	}
	normalizeMatrices(_spVec,_forceDistr);
 
 }

 Vdouble mainbb::computeFreq(codon &codonAlph){
	Vdouble pi;
	nucleotide alph;
	sequenceContainer nucSc;
	ifstream in(_options->seqfile.c_str());
	nucSc = recognizeFormat::readUnAligned(in, &alph);
	nucSc.changeGaps2MissingData();
	in.close();
	pi = freqCodonF3x4(nucSc,&codonAlph);
	makeSureNoZeroFreqs(pi);
	return pi;
}

 Vdouble mainbb::computeGTRFreq(nucleotide &nucAlph){
	Vdouble pi;
	nucleotide alph;
	sequenceContainer nucSc;
	ifstream in(_options->seqfile.c_str());
	nucSc = recognizeFormat::readUnAligned(in, &alph);
	nucSc.changeGaps2MissingData();
	in.close();
	pi = freqGTR(nucSc,&nucAlph);
	makeSureNoZeroFreqs(pi);
	return pi;
	
}
Vdouble mainbb::freqGTR(const sequenceContainer &nucSc, nucleotide * nucAlph){
	Vdouble freqGTR(nucAlph->size(),0.0);
	int pos= 0;
	int nPos = 0;
	
//	freqGTR.resize(nucSc.alphabetSize(),0.0);

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
			if ((*sIt >= 0) && (*sIt <freqGTR.size())) ++freqGTR[(*sIt)];
			if (*sIt == 4) ++freqGTR[3]; //for T (4) to U (3)
			++sIt;
			++pos;
		}
		++tIt;
	}
	changeCountsToFreqs(freqGTR);


	/*Vdouble freqGTR(nucAlph->size(),0.0);

	nucleotide n;
	for (int c = 0; c<freqGTR.size();c++){

		string s = nucAlph->fromInt(c);
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
		freqGTR[i]+=ep;
	}*/

	return freqGTR;


}

 void mainbb::replaceSequences(sequenceContainer &sc2change,sequenceContainer &originSc)
 {
	 for (int s = 0; s < originSc.numberOfSeqs();s++)
	 {
		 string name = originSc[s].name();
		 for ( int i = 0;i<sc2change.numberOfSeqs(); i++)
		 {
			 if (sc2change[i].name() == name)
			 {
				sc2change[i] = originSc[s];
				break;
			 }
		 }
			
	 }
 }
