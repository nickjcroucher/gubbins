// $Id: fromInstructionFile.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ____FROM_INSTRUCTION__FILE
#define ____FROM_INSTRUCTION__FILE

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "alphabet.h"
#include "sequenceContainer.h"
#include "someUtil.h"

#include <string>
#include <iostream>
#include <vector>
#include <map>
using namespace std;



class fromInstructionFile {
public:
	explicit fromInstructionFile(const string& instructionFileName);
	void readInstructionFile(const string& str);
	const string&searchStringInLines(const string& key) const;
	bool doesWordExistInLines(const string& key) const;
	const string& searchStringInLines(const string& key, const int index) const;
	bool getIntValueConnectedWithWord(const string& wordToSearch, int & res);
	
	
	
	void setLogFile();
	void getStartingStochasticProcess(vector<stochasticProcess>& spPtrVec,VVdouble* freqs=NULL);
	void getOneStartingStochasticProcess(stochasticProcess& sp, Vdouble * freqs = NULL);
	void getOneStartingGammaParameter(stochasticProcess& sp);
	bool getStartingEvolTrees(vector<tree>& vtree);// true if thelist tree1 file1, tree2 file2 is found.
	bool getStartingEvolTrees(vector<tree>& vtree, vector<char>& constraintsOfT0);// true if thelist tree1 file1, tree2 file2 is found.
	tree* getOneStartingEvolTree(vector<char>* constraintsOfT0);// ALOCATE NEW TREE AND NEW CONSTRAINT VECTOR.
	void getStartingSequenceData(vector<sequenceContainer>& sdPtrVec,
		const vector<alphabet* >& _alphabets);
	void getOneStartingSequenceData(sequenceContainer& sdPtrVec,
		const alphabet* _alphabets);
	void getAlphabets(vector<alphabet* >& _alphabets);// alocate with new
					// have to be deleted by the users!
	alphabet* getOneAlphabet();
	bool useGamma() {
		return doesWordExistInLines("gamma");
	}
	void getStartingGammaParameters(vector<stochasticProcess>& spPtrVec);
	void getStartingGlobalRates(vector<stochasticProcess>& spPtrVec);
	string getOutFile();
protected:

    map<string, string> _lines;
	const int _maxNumOfFiles;// = 1000;
	void getStartingGammaParameter(vector<stochasticProcess>& spPtrVec);
//	tree getStartingEvolTree();

};
#endif
