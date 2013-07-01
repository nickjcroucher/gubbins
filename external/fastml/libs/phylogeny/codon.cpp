// $Id: codon.cpp 5981 2009-03-17 14:39:39Z rubi $

#include "codon.h"
#include "nucleotide.h"
#include "amino.h"
#include "logFile.h"
#include "definitions.h"
#include "someUtil.h"
#include "matrixUtils.h"
#include "sequenceContainer.h"
#include <sstream>
#include <cctype>
#define INITIATION_CODON "i"

vector<vector<codonUtility::diffType> > codonUtility::_trtvDiff;
vector<vector<codonUtility::replacementType> > codonUtility::_synNonsynDiff;
vector<vector<codonUtility::nucDiffPlaceType> > codonUtility::_nucDiffPlace;
vector<vector<codonUtility::nucsDiffType> > codonUtility::_nucsDiff;


codon::codon(){
	geneticCodeString gcs=geneticCodeHolder::nuclearStandard;
	init(gcs);
}

codon::codon(const geneticCodeString& matrixFileString){
	init(matrixFileString);
}

void codon::init(const geneticCodeString& matrixFileString)
{
	readMatrixFromFile(matrixFileString.Val);
}

void codon::readMatrixFromFile(const string& matrixFileName){ //default value: "nuclearCode.txt"
  //	cout<<"in codon constructor"<<endl;
	stringstream in(matrixFileName.c_str());
	if (!in) {
		errorMsg::reportError("in codon::readMatrixFromFile: unable to open matrix data file");
	}

	int aa = -1; //initialized as -1 so in first iteration will change to 0
	int noOfCodons = 0;
	string strAmino; 
	bool isInitCodon = false;
	while (!in.eof()) { //20 amino acids and stop 
		string val;
		in>>val;
		if (val.size()==1) { //amino acid
			if(val == INITIATION_CODON)
				isInitCodon = true;
			else{
				aa++;
				strAmino=val;
				if (strAmino=="*") { _alphabetSize=noOfCodons;}
				isInitCodon = false;
			}
		}
		
		else if (val.size()==3 && val[0]!='#'){ //codon, # symbolizes a comment
			if(isInitCodon){
				map <string,int>::const_iterator iniItr =_codon2Int.find(val);
				if(iniItr == _codon2Int.end())
					errorMsg::reportError("Initiation codon with undefined index at codon::readMatrixFromFile");
				else
					_initiationIndex2codon[iniItr->second] = val;
			}
			else{
				_geneticCode[val]=strAmino;
				_codon2Int[val]=noOfCodons;
				noOfCodons++;
			}
		}
		else {
			
			if (noOfCodons!=64){
				string err="in codon::readMatrixFromFile: total number of codons = "+int2string(noOfCodons);
				errorMsg::reportError(err);
			}
			return;
		}
	}
}
codon& codon::operator=(const codon& other) {
	_geneticCode = other._geneticCode; //key - codon, value - amino acid
	_codon2Int = other._codon2Int;//key string of codon int= integer value of codon
	_alphabetSize = other._alphabetSize;
	_initiationIndex2codon = other._initiationIndex2codon;
	return *this;
}
// codon::codon(const  codon& other):
// 	_geneticCode(other._geneticCode), //key - codon, value - amino acid
// 	_codon2Int(other._codon2Int),//key string of codon int= integer value of codon
// 	_alphabetSize(other._alphabetSize){}


//return -99 if not succeeds.
int codon::fromChar(const string& s, const int pos) const {
	if (s.size() <= pos+2) {
		//errorMsg::reportError("Trying to read a codon pass the end of the string. The number of nucleotide may not be divisible by three");
		string textToPrint("Trying to read a codon pass the end of the string. The number of nucleotide may not be divisible by three");
		LOG(1,<<textToPrint<<endl);
		return -99;
	}

	nucleotide nuc;
	int p1,p2,p3;
	p1 = nuc.fromChar(s[pos]);
	p2 = nuc.fromChar(s[pos+1]);
	p3 = nuc.fromChar(s[pos+2]);


	if ((p1 <0) || (p2 <0) || (p3 <0)) 
		return gap(); 
	else if ((p1 ==15) || (p2 ==15) || (p3 ==15)) return unknown(); // unknown.
	else if ((p1 >4) || (p2 >4) || (p3 >4)) return unknown(); //unknown.
	string strCodon="";
	//change U --> T
	if (p1==4) strCodon+="T";
	else  strCodon+=toupper(s[pos]);
	if (p2==4) strCodon+="T";
	else  strCodon+=toupper(s[pos+1]);
	if (p3==4) strCodon+="T";
	else  strCodon+=toupper(s[pos+2]);
	//const string strCodon = s.substr(pos,3);
	map <string,int> tmpMap=_codon2Int;
	map <string,int>::iterator it1;
	it1=tmpMap.find(strCodon);
	if (it1==tmpMap.end()){
		
		string err="error in codon::fromChar cannot find codon "+strCodon;
		errorMsg::reportError(err);
	}
	return tmpMap[strCodon];
}

vector<int> codon::fromString(const string &str) const {
	vector<int> vec;
	if (str.size()%3!=0) {
		errorMsg::reportError("error in function codon::fromString. String length should be a multiplication of 3");
	}
	for (int i=0;i<str.size();i+=3)
	  vec.push_back(fromChar(str,i));
	return vec;
}

string codon::fromInt(const int in_id)  const{
	if (in_id == unknown())
		return "XXX";
	if (in_id == gap()) 
		return "---";
	map <string, int> tmpMap = _codon2Int;
	map <string, int>::iterator it=tmpMap.begin();
	while (it!=tmpMap.end()){
		if ((*it).second==in_id){
			return (*it).first;
		}
		it++;
	}
	string err="error in function codon::fromInt: no codon found for the integer";
	errorMsg::reportError(err);
	return (string("we should never get here - the reportError above will exit"));
}

codonUtility::replacementType codonUtility::codonReplacement(const int c1, const int c2, const codon &cod){
	if (c1 == c2) return codonUtility::sameCodon;
	else if (codonUtility::aaOf(c1,cod) == codonUtility::aaOf(c2,cod)) return codonUtility::synonymous;
	return codonUtility::non_synonymous;
}

int codonUtility::aaOf(const int c1, const codon &cod){
    amino a;
    if (c1==cod.gap()) 
		return a.gap();
	if (c1==cod.unknown())
		return a.unknown();
	string strCodon=cod.fromInt(c1);
	map <string,string> geneticCode=cod.geneticCode();
	map <string,string>::iterator pos;
	if ((pos=geneticCode.find(strCodon)) == geneticCode.end()){
		string err="error in codonUtility::aaOf: cannot find codon "+strCodon;
		errorMsg::reportError(err);
	}
	if (pos->second.size() > 1){
		errorMsg::reportError("error in codonUtility::aaOf: amino acid 1 letter code > 1");
	}
	return a.fromChar(*pos->second.c_str());
}


codonUtility::diffType codonUtility::codonDiff(const int c1, const int c2, const codon &cod){
	if (c1==c2) return codonUtility::equal;
	nucleotide n;
	string s1 = cod.fromInt(c1);
	string s2 = cod.fromInt(c2);

	int pos1 = n.fromChar(s1[0])+n.fromChar(s2[0]);
	int pos2 = n.fromChar(s1[1])+n.fromChar(s2[1]);
	int pos3 = n.fromChar(s1[2])+n.fromChar(s2[2]);

	if (s1[0]!=s2[0] && s1[1]!=s2[1] && s1[2]!=s2[2])
		return  codonUtility::threesub; 

	if (s1[0]==s2[0] && s1[1]==s2[1] && s1[2]!=s2[2]) {
		if (pos3%2==0) return codonUtility::tr;
		else return codonUtility::tv;
	}
	if (s1[1]==s2[1] && s1[2]==s2[2] && s1[0]!=s2[0]) {
		if (pos1%2==0) return codonUtility::tr;
		else return codonUtility::tv;
	}
	if (s1[0]==s2[0] && s1[2]==s2[2] && s1[1]!=s2[1]) {
		if (pos2%2==0) return codonUtility::tr;
		else return codonUtility::tv;
	}

	if (s1[0]==s2[0] && pos2%2==0 && pos3%2==0)
		return  codonUtility::twoTrs;
	if (s1[1]==s2[1] && pos1%2==0 && pos3%2==0)
		return  codonUtility::twoTrs;
	if (s1[2]==s2[2] && pos1%2==0 && pos2%2==0)
		return  codonUtility::twoTrs;

	if (s1[0]==s2[0] && pos2%2!=0 && pos3%2!=0)
		return  codonUtility::twoTvs;
	if (s1[1]==s2[1] && pos1%2!=0 && pos3%2!=0)
		return  codonUtility::twoTvs;
	if (s1[2]==s2[2] && pos1%2!=0 && pos2%2!=0)
		return  codonUtility::twoTvs;

	return codonUtility::trtv;
}


//return the place (0, 1, or 2) that the two codons are different 
//and the identity of the different nucleotide in the target codon.
//For example, nucDiffPlace(ATG, ACG) retruns C2
codonUtility::nucDiffPlaceType codonUtility::nucDiffPlace(const int fromCodon, const int targetCodon, const codon &cod){
	if (fromCodon == targetCodon) 
		return codonUtility::EQUAL;

	codonUtility::nucDiffPlaceType res = A1;
	nucleotide nuc;
	string s1 = cod.fromInt(fromCodon);
	string s2 = cod.fromInt(targetCodon);

	int diffNum = 0;
	if (s1[0] != s2[0]){
		++diffNum;
		switch (s2[0])
		{
		case 'A': res = A1;
			break;
		case 'C': res = C1;
			break;
		case 'G': res = G1;
			break;
		case 'T': res = T1;
			break;
		default: 
			errorMsg::reportError("error in codonUtility::nucDiffPlace.");
			break;
		}
	}
	if (s1[1] != s2[1]){
		++diffNum;
		switch (s2[1])
		{
		case 'A': res = A2;
			break;
		case 'C': res = C2;
			break;
		case 'G': res = G2;
			break;
		case 'T': res = T2;
			break;
		default: 
			errorMsg::reportError("error in codonUtility::nucDiffPlace.");
			break;
		}
	}
	if (s1[2] != s2[2]){
		++diffNum;
		switch (s2[2])
		{
		case 'A': res = A3;
			break;
		case 'C': res = C3;
			break;
		case 'G': res = G3;
			break;
		case 'T': res = T3;
			break;
		default: 
			errorMsg::reportError("error in codonUtility::nucDiffPlace.");
			break;
		}
	}
	if (diffNum == 0)
		errorMsg::reportError("error in codonUtility::nucDiffPlace. Can't find different nucleotide");
	if (diffNum > 1)
		res = MUL_SUB;
	return res;
}

//return the different nucleotides between the fron and target codons.
//For example, nucsPlace(ATG, ACG) retruns TC
codonUtility::nucsDiffType codonUtility::nucsDiff(const int fromCodon, const int targetCodon, const codon &cod){
	if (fromCodon == targetCodon) 
		return codonUtility::SAME;

	codonUtility::nucsDiffType res = AC;
	nucleotide nuc;
	string s1 = cod.fromInt(fromCodon);
	string s2 = cod.fromInt(targetCodon);

	int diffNum = 0;
	int from = 0;
	int to = 0;
	if (s1[0] != s2[0])
	{	
		++diffNum;
		from = s1[0];
		to = s2[0];
	}
	if (s1[1] != s2[1])
	{	
		++diffNum;
		from = s1[1];
		to = s2[1];
	}
	if (s1[2] != s2[2])
	{
		++diffNum;
		from = s1[2];
		to = s2[2];
	}
	switch(from)
	{
		case 'A':
			switch(to)
			{
			case 'G':res = AG;break;
			case 'T':res = AT;break;
			case 'C':res = AC;break;
			default: 
				errorMsg::reportError("error in codonUtility::nucsDiff.");
				break;
			}
			break;
		case 'G':
			switch(to)
			{
			case 'A':res = AG;break; 
			case 'T':res = GT;break;
			case 'C':res = CG;break;
			default: 
				errorMsg::reportError("error in codonUtility::nucsDiff.");
				break;
			}
			break;
		case 'C':
			switch(to)
			{
			case 'G':res = CG;break; 
			case 'T':res = CT;break;
			case 'A':res = AC;break;
			default: 
				errorMsg::reportError("error in codonUtility::nucsDiff.");
				break;
			}
			break;
		case 'T':
			switch(to)
			{
			case 'G':res = GT;break; 
			case 'A':res = AT;break;
			case 'C':res = CT;break;
			default: 
				errorMsg::reportError("error in codonUtility::nucsDiff.");
				break;
			}
			break;
		default: 
			errorMsg::reportError("error in codonUtility::nucsDiff.");
			break;
	}
	
	if (diffNum == 0)
		errorMsg::reportError("error in codonUtility::nucsDiff. Can't find different nucleotide");
	if (diffNum > 1)
		res = DIFF;
	return res;
}



void codonUtility::initSubMatrices(const codon& cod){

	if ((_trtvDiff.size() == cod.size()) && (_synNonsynDiff.size() == cod.size()) && (_nucDiffPlace.size() == cod.size()) && (_nucsDiff.size() == cod.size()))
		return;

	_trtvDiff.resize(cod.size());
	_synNonsynDiff.resize(cod.size());
	_nucDiffPlace.resize(cod.size());
	_nucsDiff.resize(cod.size());
	for (int i = 0; i < _trtvDiff.size(); ++i)
	{
		_trtvDiff[i].resize(cod.size());
		_synNonsynDiff[i].resize(cod.size());
		_nucDiffPlace[i].resize(cod.size());
		_nucsDiff[i].resize(cod.size());

	}
	//resizeMatrix<diffType>(_trtvDiff, cod.size(), cod.size());
	//resizeMatrix<replacementType>(_synNonsynDiff, cod.size(), cod.size());
	//resizeMatrix<nucDiffPlaceType>(_nucDiffPlace, cod.size(), cod.size());
	for (int i = 0; i < cod.size(); ++i){
		for (int j =0; j <= i; ++j){
			_trtvDiff[i][j] = _trtvDiff[j][i] = codonDiff(i, j, cod);
			_synNonsynDiff[i][j] = _synNonsynDiff[j][i] = codonReplacement(i, j, cod);
			_nucDiffPlace[i][j] = nucDiffPlace(i, j, cod);
			_nucDiffPlace[j][i] = nucDiffPlace(j, i, cod);
			_nucsDiff[i][j] =  nucsDiff(i,j,cod);
			_nucsDiff[j][i] = nucsDiff(j,i,cod);
		}
	}
}

//returns the number (codonCounter) and frequency (codonUsage) of each codon in the sequnece container
void codonUtility::getCodonUsage(const sequenceContainer& sc, Vint& codonCounter, Vdouble& codonUsage)
{
	if (sc.getAlphabet()->size() != 61)
		errorMsg::reportError("cannot calculate codon usage when alphabet is not codon");
	codonCounter.resize(61, 0);
	codonUsage.resize(61, 0.0);
	codon alph;
	int sum = 0;
	for (int s = 0; s < sc.numberOfSeqs();++s) {
		int id = sc.placeToId(s);
		for (int pos = 0; pos < sc.seqLen(); ++pos)
		{
			int cod = sc[id][pos];
			if (alph.isSpecific(cod))
			{
				++sum;
				++codonCounter[cod];
			}
		}
	}

	for (int c = 0; c < codonCounter.size(); ++c)
		codonUsage[c] = static_cast<MDOUBLE>(codonCounter[c]) / sum;
}


//in codonUsageFile: only 3-letter-codon and frequency seperated by "\t"
void codonUtility::readCodonUsage(const string& codonUsageFileName, Vdouble& codonUsage,const codon &alph)
{
	codonUsage.resize(alph.size(), 0.0);
	ifstream inFile(codonUsageFileName.c_str());
	vector<string> inFileData;
	putFileIntoVectorStringArray(inFile, inFileData);
	inFile.close();
	if (inFileData.empty()){
		errorMsg::reportError("unable to open file, or file is empty in codonUtility::readCodonUsage");
	}

	vector<string>::const_iterator it = inFileData.begin();
	for (; it!= inFileData.end(); ++it) 
	{
		if (it->empty()) //empty line
			continue; 
		int endCodon = it->find_first_of("\t", 0);
		int startFreq = it->find_first_not_of("\t ", endCodon);
        if (startFreq>0)
		{
			string codonStr = it->substr(0, endCodon);
			string freqStr = it->substr(startFreq);
			MDOUBLE freq = string2double(freqStr);
			if(freq == 0.0) freq = EPSILON;
			codonUsage[alph.fromChar(codonStr, 0)] = freq;
		}
	}
}

//calculates the CAI for the whole MSA and for each position.
//The calculation is based on a pre-calculated codonUsage vector.
//The calculation is based on Sharp & Li (1987) NAR, 15:1281-1295
MDOUBLE codonUtility::calcCodonAdaptationIndex(const sequenceContainer& sc, const Vdouble& codonUsage, Vdouble& cai4site)
{	
	//the returned value: calculated as the average CAI for the MSA, rather than the geometrical mean as in Sharp & Li
	MDOUBLE wholeAlignmentCai = 0.0; 
	codon alph;
	amino am;
	//1. calculate Wk = the frequency of codon k relative to the frequency of the optimal codon for that amino acid.
	Vdouble Wk(codonUsage.size(), 0.0);
	int aaId;
	for (aaId = 0; aaId < am.size(); ++aaId)
	{
		Vint codonsOfAa = aminoUtility::codonOf(aaId, alph);
		//finding the most frequent codon for this aa
		MDOUBLE mostFrequent = 0.0;
		Vint::const_iterator iter;
		for (iter = codonsOfAa.begin(); iter != codonsOfAa.end(); ++iter)
		{
			if (codonUsage[*iter] > mostFrequent)
				mostFrequent = codonUsage[*iter];
		}

		//calculating Wk
		for (iter = codonsOfAa.begin(); iter != codonsOfAa.end(); ++iter)
			Wk[*iter] = codonUsage[*iter] / mostFrequent;
	}

	//2. calculate CAI
	cai4site.resize(sc.seqLen(), 0.0);
	int pos;
	for (pos = 0; pos < sc.seqLen(); ++pos)
	{
		MDOUBLE cai = 0.0;
		int informativeCodons = 0;
		for (int s = 0; s < sc.numberOfSeqs();++s) 
		{
			int id = sc.placeToId(s);
			int cod = sc[id][pos];
			if(!alph.isSpecific(cod))
				continue;
			cai += Wk[cod];
			++informativeCodons;
		}
		
		cai /= static_cast<MDOUBLE>(informativeCodons);
		cai4site[pos] = cai;
		wholeAlignmentCai += cai;
	}
	return wholeAlignmentCai;
}



bool codon::isStopCodon(const int in_id) const
{
	if (in_id == unknown()) return false;
	if (in_id == gap()) return false;
	if ((in_id >= 0 ) && (in_id < _alphabetSize)) return false;
	return true;
}

bool codon::isInitiationCodon(const int in_id) const
{
	bool result = true;
	map <int,string>::const_iterator itr = _initiationIndex2codon.find(in_id);
	if(itr == _initiationIndex2codon.end()){
		result = false;
	}
	return result;
}


