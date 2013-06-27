// $Id: fromInstructionFile.cpp 962 2006-11-07 15:13:34Z privmane $

#include "definitions.h"
#include "fromInstructionFile.h"
#include "treeUtil.h"
#include "nucleotide.h"
#include "amino.h"
#include "uniDistribution.h"
#include "gammaDistribution.h"
#include "readDatMatrix.h"
#include "aaJC.h"
#include "nucJC.h"
#include "hky.h"
#include "trivialAccelerator.h"
#include "chebyshevAccelerator.h"
#include "phylipFormat.h"
#include "maseFormat.h"
#include "fastaFormat.h"
#include "clustalFormat.h"
#include "molphyFormat.h"
#include "datMatrixHolder.h"
#include "someUtil.h"

#include <iostream>
#include <fstream>
#include <memory>
#include <iterator>
#include <cstdio>
using namespace std;

//#define VERBOS

void fromInstructionFile::readInstructionFile(const string& str){
	ifstream f;
	f.open(str.c_str());
	if (f==NULL) {
	  string tmp = "Unable to open the instraction file : \""+str+"\""; 
	  errorMsg::reportError(tmp);
	}
	string key, value;
	while (!f.eof()){
		f >> key;
		if (!key.empty()){
			toLower(key);// put the key in lower case.
			getline(f,value);
			value.erase(0,value.find_first_not_of(" \t")); // clear leading white space 
			_lines[key]=value;
		}
	}
	f.close();
}

fromInstructionFile::fromInstructionFile(const string& str):_maxNumOfFiles(1000){
  readInstructionFile(str);
}

// THIS IS NOT WORKING ON SOME OLD VERSIONS OF g++
//string I2A(const int & v)
//{
//  stringstream s("");
//  s<<v;
//  return(s.str());
//}
//
//string F2A(const float & v)
//{
//  stringstream s("");
//  s<<v;
//  return(s.str());
//}

string I2A(const int & v)
{
  char buf[100];
  sprintf(buf,"%d",v);
  return buf;
}

string F2A(const float & v)
{
  char buf[100];
  sprintf(buf,"%f",v);
  return buf;
}




bool fromInstructionFile::doesWordExistInLines(const string& key) const{
  return (_lines.count(key)>0);
}

const string & fromInstructionFile::searchStringInLines(const string& key) const
{
#ifdef VERBOS
	map<string, string>::const_iterator pos;
	pos = _lines.begin();
	for (; pos != _lines.end(); ++pos) {
        cout << "key: \"" << pos->first << "\" "
             << "value: " << pos->second << endl;
    }
#endif
	
	
	
  static const string emptystr("");
  if (_lines.count(key) > 0)
    return(_lines.find(key)->second);
  else
    return(emptystr);
}

const string& fromInstructionFile::searchStringInLines(const string& key, const int index) const
{
  static const string emptystr("");

  string realKey(key+int2string(index));

  if (_lines.count(realKey) > 0)
    return(_lines.find(realKey)->second);
  else
    return(emptystr);
}

void fromInstructionFile::setLogFile() {
	string logfilename(searchStringInLines("logfile"));
	if (logfilename == "") logfilename = "-";

	if (logfilename == "-") {
	  myLog::setLogOstream(&cout);
	}
	else{
	  ofstream* outLF = new ofstream(logfilename.c_str());
	  if (!outLF) {
	    errorMsg::reportError("unable to open file for reading");
	  }
	  myLog::setLogOstream(outLF);
	}
	string loglvl(searchStringInLines("loglvl"));
	if (loglvl=="") myLog::setLogLvl(3); // default value
	else myLog::setLogLvl(atoi(loglvl.c_str()));
	LOG(3,<<"START OF LOG FILE\n\n");
}

bool fromInstructionFile::getIntValueConnectedWithWord(const string& wordToSearch,
													   int & val){
  string p(searchStringInLines(wordToSearch));	
	if (p == "") {
		return false;
	}
	val=atoi(p.c_str());
	return true;
}

string fromInstructionFile::getOutFile() {
	string outfilename(searchStringInLines("outfile"));
	if (outfilename == "") outfilename = "-";
	return outfilename;
}

void fromInstructionFile::getAlphabets(vector<alphabet* >& _alphabets) {
	if (_alphabets.size() !=0) {errorMsg::reportError("error in fromInstructionFile::getAlphabetSize");}
	for (int i=1; i < _maxNumOfFiles; ++i ) {
		string p(searchStringInLines("alphabet",i));
		if (p == "") return;
		int alphRes = atoi(p.c_str());
		if (alphRes == 4) {
			alphabet* alp = new nucleotide;
			_alphabets.push_back(alp);
		}
		else if (alphRes == 20) {
			alphabet* alp = new amino;
			_alphabets.push_back(alp);
		}
		else errorMsg::reportError("No relaven number after the word alphabet in the instruction file.");
	}
	for (size_t z=1; z< _alphabets.size(); ++z) {
		if (_alphabets[z]!= _alphabets[0]) {
			errorMsg::reportError("currently all seq. must be of the same alphabet size");
		}
	}
}

alphabet* fromInstructionFile::getOneAlphabet( ) {
	alphabet* _alphabet = NULL;
	int alphRes;
	
	bool ok = getIntValueConnectedWithWord("alphabet",alphRes);
	if (!ok) {
		ok = getIntValueConnectedWithWord("alphabet1",alphRes);
	
		if (!ok) errorMsg::reportError("didn't find alphabet size in instruction file");
	}if (ok==true) {
		if (alphRes == 4) {
			_alphabet = new nucleotide;
		}
		else if (alphRes == 20) {
			_alphabet = new amino;
		}
		else errorMsg::reportError("No number after the word alphabet in the instruction file.");
	}
	return _alphabet;
}

void fromInstructionFile::getOneStartingStochasticProcess(stochasticProcess& sp, Vdouble * freqs){
	bool useGamma = doesWordExistInLines("gamma");
	distribution *dist = NULL;
	if (!useGamma) dist =  new uniDistribution;
	else dist =  new gammaDistribution(1,4);
		
	replacementModel *probMod=NULL;
	pijAccelerator *pijAcc=NULL;

	string wordUse = "model";
	bool usemodel1 = doesWordExistInLines("model1");
	if (usemodel1 == true) wordUse="model1";

	string modelName(searchStringInLines(wordUse));// we can use model or model1
	if (modelName == "") {
		errorMsg::reportError("could not find model name in instruction file");
	}
	
	if (strcmp(modelName.c_str(),"day")==0) {
		(freqs==NULL)? probMod=new pupAll(datMatrixHolder::dayhoff) : probMod=new pupAll(datMatrixHolder::dayhoff,*freqs);
		pijAcc = new chebyshevAccelerator(probMod); 
	}
	else if (strcmp(modelName.c_str(),"jtt")==0) {
		(freqs==NULL)? probMod=new pupAll(datMatrixHolder::jones):probMod=new pupAll(datMatrixHolder::jones,*freqs) ;
		pijAcc =new chebyshevAccelerator(probMod);
	}
	else if (strcmp(modelName.c_str(),"rev")==0) {
		(freqs==NULL)? probMod=new pupAll(datMatrixHolder::mtREV24) : probMod=new pupAll(datMatrixHolder::mtREV24,*freqs);
		pijAcc = new chebyshevAccelerator(probMod);
	}
	else if (strcmp(modelName.c_str(),"wag")==0) {
		(freqs==NULL)? probMod=new pupAll(datMatrixHolder::wag) : probMod=new pupAll(datMatrixHolder::wag, *freqs);
		pijAcc = new chebyshevAccelerator(probMod);
	}
	else if (strcmp(modelName.c_str(),"cprev")==0) {
		(freqs==NULL)?  probMod=new pupAll(datMatrixHolder::cpREV45) : probMod=new pupAll(datMatrixHolder::cpREV45, *freqs);
		pijAcc = new chebyshevAccelerator(probMod);
	}
	else if (strcmp(modelName.c_str(),"nucjc")==0) {
		probMod=new nucJC; pijAcc = new trivialAccelerator(probMod);
	}
	else if (strcmp(modelName.c_str(),"aaJC")==0) {
		probMod=new aaJC; pijAcc = new trivialAccelerator(probMod);
	}
	else if (modelName=="hky"||modelName=="k2p") {
	  MDOUBLE ratio (atof(searchStringInLines("ratio").c_str())); // get alpha
	  MDOUBLE Ap(0.25), Cp(0.25), Gp(0.25), Tp(0.25);
	  sscanf(searchStringInLines("ACGprob").c_str(),"%lf,%lf,%lf", &Ap, &Cp, &Gp);
	  Tp=1.0-(Ap+Cp+Gp);
	  probMod=new hky(Ap,Cp,Gp,Tp,ratio); pijAcc = new trivialAccelerator(probMod);
	}
	else {
		errorMsg::reportError("This replacement model is not yet available");
	}

	stochasticProcess s1s(dist, pijAcc);
	if (probMod) delete probMod;
	if (pijAcc) delete pijAcc;
	if (dist) delete dist;
	sp = s1s;
}

void fromInstructionFile::getStartingStochasticProcess(vector<stochasticProcess>& spPtrVec, VVdouble* freqs) {
	if (spPtrVec.size() !=0) {errorMsg::reportError("error in fromInstructionFile::getStartingSequenceData");}
	bool useGamma = doesWordExistInLines("gamma");
	for (int i=0; i < _maxNumOfFiles; ++i) {
		Vdouble* freq_i = (freqs==NULL) ? NULL: &((*freqs)[i]);

		distribution *dist = NULL;
		if (!useGamma) dist =  new uniDistribution;
		else dist =  new gammaDistribution(1,4);
		

		replacementModel *probMod=NULL;
		pijAccelerator *pijAcc=NULL;
		string model(searchStringInLines("model",i+1));
		if (model == "") return;
		if (model=="day") {
			if (freq_i == NULL) {
				probMod=new pupAll(datMatrixHolder::dayhoff);//pijAcc = new chebyshevAccelerator(probMod); 
			} else {
				probMod=new pupAll(datMatrixHolder::dayhoff,*freq_i);//pijAcc = new chebyshevAccelerator(probMod); 
			}
			pijAcc = new trivialAccelerator(probMod);
		}
		else if (model=="jtt") {
			if (freq_i == NULL) {
				probMod=new pupAll(datMatrixHolder::jones) ; //pijAcc =new chebyshevAccelerator(probMod);
			}
			else {
					probMod=new pupAll(datMatrixHolder::jones,*freq_i) ; //pijAcc =new chebyshevAccelerator(probMod);
			}
			pijAcc = new trivialAccelerator(probMod);
		}
		else if (model=="rev") {
			if (freq_i == NULL) {
				probMod=new pupAll(datMatrixHolder::mtREV24);//pijAcc = new chebyshevAccelerator(probMod);
			} else {
					probMod=new pupAll(datMatrixHolder::mtREV24,*freq_i);//pijAcc = new chebyshevAccelerator(probMod);
			}
			pijAcc = new trivialAccelerator(probMod);
		} else if (model=="wag") {
			if (freq_i == NULL) {
				probMod=new pupAll(datMatrixHolder::wag);//pijAcc = new chebyshevAccelerator(probMod);
			} else {
				probMod=new pupAll(datMatrixHolder::wag,*freq_i);//pijAcc = new chebyshevAccelerator(probMod);
			}
			pijAcc = new trivialAccelerator(probMod);
		}	else if (model=="cprev") {
			if (freq_i == NULL) {
				probMod=new pupAll(datMatrixHolder::cpREV45);//pijAcc = new chebyshevAccelerator(probMod);
			} else {
				probMod=new pupAll(datMatrixHolder::cpREV45,*freq_i);//pijAcc = new chebyshevAccelerator(probMod);
			}
			pijAcc = new trivialAccelerator(probMod);
		}
		else if (model == "nucjc") {
			probMod=new nucJC; pijAcc = new trivialAccelerator(probMod);
		}
		else if (model == "aaJC") {
			probMod=new aaJC; pijAcc = new trivialAccelerator(probMod);
		}
		else {errorMsg::reportError("This replacement model is not yet available");
		}

		stochasticProcess s1s(dist, pijAcc);
		spPtrVec.push_back(s1s);
		if (probMod) delete probMod;
		if (pijAcc) delete pijAcc;
		if (dist) delete dist;
	}
}

bool fromInstructionFile::getStartingEvolTrees(vector<tree>& vtree,vector<char>& constraintsOfT0){
	if (vtree.size() !=0) {
		errorMsg::reportError("error in fromInstructionFile::getStartingEvolTrees");
	}
	string oneTreeFileName(searchStringInLines("treefile"));
	if (oneTreeFileName =="" ) {
		errorMsg::reportError("The tree file name must be given in the instruction file");
	}
	getStartingTreeVecFromFile(oneTreeFileName,vtree,constraintsOfT0);
	for (size_t k=0;k<vtree.size();++k) {
		if (!vtree[k].withBranchLength()) vtree[k].createFlatLengthMatrix(0.05);
	}
	return true;
}


bool fromInstructionFile::getStartingEvolTrees(vector<tree>& vtree){
	if (vtree.size() !=0) {errorMsg::reportError("error in fromInstructionFile::getStartingEvolTrees");}
//	for (int i=1; i < _maxNumOfFiles; ++i ) {
//		auto_ptr<string> treeFileName(searchStringInFile("treefile",i,_instructionFile));
//		if ((treeFileName.get() == NULL) && (i==1)) {
			string oneTreeFileName(searchStringInLines("treefile"));
			if (oneTreeFileName=="" ) {
				errorMsg::reportError("The tree file name must be given in the instruction file");
			}
			vtree = getStartingTreeVecFromFile(oneTreeFileName);
				//tree tmpT(*oneTreeFileName);
				//vtree.push_back(tmpT);
			for (size_t k=0;k<vtree.size();++k) {
				if (!vtree[k].withBranchLength()) 
					vtree[k].createFlatLengthMatrix(0.05);
			}
			return true;
//		}
//		if (treeFileName.get() == NULL) return true;// found some trees
//		tree t1(*treeFileName);
//		if (!t1.WithBranchLength()) t1.create_flat_length_matrix(0.05);
//		vtree.push_back(t1);
//	}
//	errorMsg::reportError("error in function fromInstructionFile::getStartingEvolTrees");
//	return false;
}

void fromInstructionFile::getStartingSequenceData(vector<sequenceContainer>& sdPtrVec,
												  const vector<alphabet* >& _alphabets){
	if (sdPtrVec.size() !=0) {errorMsg::reportError("error in fromInstructionFile::getStartingSequenceData");}
	for (int i=1; i <= _maxNumOfFiles; ++i ) {
		string sequenceFileName(searchStringInLines("seqfile",i));
		if ((sequenceFileName == "") && (i==1)) sequenceFileName="-";
		else if (sequenceFileName == "") return;

		istream* inPtr;
		if (sequenceFileName == "-") {
		  LOG(5,<<"in this option, the sequences are inputed from cin\n...");
		  inPtr = &cin;
		}else{
		  inPtr = new ifstream(sequenceFileName.c_str());
		}
		istream& in = *inPtr;
		sequenceContainer original;

		string sequenceFileFormat(searchStringInLines("format",i));
		if ((sequenceFileFormat == "") && (i>1)) {// it is probably the format of number 1.
			string sequenceFileFormatOf1(searchStringInLines("format",1));
			sequenceFileFormat = sequenceFileFormatOf1;
		}
		alphabet* currentAlphabet = NULL;
		if ((_alphabets.size() == 1) && (i > 1)) currentAlphabet = _alphabets[0];
		else {
			currentAlphabet = _alphabets[i-1];
		}
		if      (sequenceFileFormat== "mase")	original= maseFormat::   read(in,currentAlphabet);
		else if (sequenceFileFormat=="molphy")  original= molphyFormat:: read(in,currentAlphabet);
		else if (sequenceFileFormat=="clustal") original= clustalFormat::read(in,currentAlphabet);
		else if (sequenceFileFormat=="fasta")	original= fastaFormat::  read(in,currentAlphabet);
		else if (sequenceFileFormat=="phylip")	original= phylipFormat:: read(in,currentAlphabet);
		else errorMsg::reportError(" format not implemented yet in this version... ");
		
//		if (original == NULL) errorMsg::reportError(" unable to find/open input sequence file");
		
		if (doesWordExistInLines("removeGapPositions")) {
//			vector<int> parCol;
//			original.getParticiantColVecAccordingToGapCols(parCol);
//			sequenceData _sd(*original,parCol);
//			sdPtrVec.push_back(_sd);
//			delete original;
			errorMsg::reportError("remove gap position is not implemented yet");
		} //else if (doesWordExistInLines("gapsToMissingData")) {
			//LOG(5,<<"gaps are changed to missing data..."<<endl);
			original.changeGaps2MissingData();
			sdPtrVec.push_back(original);
		//}
	}

}

tree* fromInstructionFile::getOneStartingEvolTree(vector<char>* constraintsOfT0) {
	tree* _tree = NULL;
	
	string wordUse = "treefile";
	bool usetreefile1 = doesWordExistInLines("treefile1");
	if (usetreefile1 == true) wordUse="treefile1";

	string treeFileName(searchStringInLines(wordUse)); // either treefile or treefile1 is OK.
	if (treeFileName=="" ) {
	  _tree = NULL;
	  constraintsOfT0 = NULL;
	  return _tree;
	}

	vector<char> constraints;
	_tree = new tree(treeFileName,constraints);
	constraintsOfT0 = new vector<char>(constraints);
	return _tree;
}

void fromInstructionFile::getOneStartingSequenceData(sequenceContainer& sd,
													 const alphabet* _alphabets) {
	ifstream ins;
	istream* inPtr = NULL;	

	string wordUse = "seqfile";
	bool useseqfile1 = doesWordExistInLines("seqfile1");
	if (useseqfile1 == true) wordUse="seqfile1";

	string sequenceFileName(searchStringInLines(wordUse)); // so it can be used with both seqfile and seqfile1
	if (sequenceFileName == "") sequenceFileName="-";
	if (sequenceFileName == "-") {
	  inPtr = &cin;
	}
	else{
	  ins.open(sequenceFileName.c_str());
	  if (! ins.is_open())
	    errorMsg::reportError("can not open sequace file");
	  inPtr = &ins;
	}
	
	istream& in = *inPtr;
	sequenceContainer original;
	
	wordUse = "format";
	bool useFormat1 = doesWordExistInLines("format1");
	if (useFormat1 == true) wordUse="format1";

	string sequenceFileFormat(searchStringInLines(wordUse));
	if (sequenceFileFormat == "") {
	  sequenceFileFormat = "fasta"; // default
	}

	if      (sequenceFileFormat == "mase")    original= maseFormat::read(in,_alphabets);
	else if (sequenceFileFormat == "molphy")  original= molphyFormat::read(in,_alphabets);
	else if (sequenceFileFormat == "clustal") original= clustalFormat::read(in,_alphabets);
	else if (sequenceFileFormat == "fasta")   original= fastaFormat::read(in,_alphabets);
	else if (sequenceFileFormat == "phylip")  original= phylipFormat::read(in,_alphabets);
	else errorMsg::reportError(" format not implemented yet in this version... ");
		
	if (doesWordExistInLines("removeGapPositions")) {
	  errorMsg::reportError("remove gap position is not implemented yet");
	} 
	//LOG(5,<<"gaps are changed to missing data..."<<endl);
	original.changeGaps2MissingData();
	sd = original;
}

void fromInstructionFile::getStartingGammaParameters(vector<stochasticProcess>& spPtrVec) {
	for (size_t i=0; i < spPtrVec.size(); ++i) {
		string alphaParam(searchStringInLines("alpha",i+1));
		if ((alphaParam == "") && (i==0)) {
		  getStartingGammaParameter(spPtrVec);
		  return;
		}
		if (alphaParam == "") {
		  MDOUBLE alpha = atof(alphaParam.c_str());
		  (static_cast<gammaDistribution*>(spPtrVec[i].distr()))->setAlpha(alpha);
		}
	}
}

void fromInstructionFile::getOneStartingGammaParameter(stochasticProcess& sp) {
	MDOUBLE alpha = 0;
	string alphaParam0(searchStringInLines("alpha",0));
	if (alphaParam0 != "") {
	  alpha = atof(alphaParam0.c_str());
	} else {
	  string alphaParam1(searchStringInLines("alpha",1));
	  if (alphaParam1 != "") {
	    alpha = atof(alphaParam1.c_str());
	  } else {
	    string alphaParam2(searchStringInLines("alpha"));
	    if (alphaParam2 != "") {
	      alpha = atof(alphaParam2.c_str());
	    } else { // no alpha parameter given,
	      return;
	    }
	  }
	}
	(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(alpha);
}

void fromInstructionFile::getStartingGammaParameter(vector<stochasticProcess>& spPtrVec) {
	string alphaParam(searchStringInLines("alpha"));
	for (size_t i=0; i < spPtrVec.size(); ++i) {
	  if (alphaParam != "") {
	    MDOUBLE alpha = atof(alphaParam.c_str());
	    (static_cast<gammaDistribution*>(spPtrVec[i].distr()))->setAlpha(alpha);
		}
	}
}

void fromInstructionFile::getStartingGlobalRates(vector<stochasticProcess>& spPtrVec) {
	for (size_t i=0; i < spPtrVec.size(); ++i) {
		string rate(searchStringInLines("rate",i+1));
		if (rate != "") {
		  MDOUBLE grate = atof(rate.c_str());
		  spPtrVec[i].setGlobalRate(grate);
		}
	}
}
