#ifndef ___SIMULATE_CODONS_JUMPS__
#define ___SIMULATE_CODONS_JUMPS__

#include "simulateJumpsAbstract.h"
using namespace std;

/******************************************************************
This class implements simulateJumpsAbstract for small alphabets: (tested so far up to 3)
*******************************************************************/

class simulateCodonsJumps:public simulateJumpsAbstract  {
public:
	simulateCodonsJumps(const tree& inTree, const stochasticProcess& sp, const int alphabetSize);
	virtual ~simulateCodonsJumps();
	
	//for a branch length specified by a nodeName: 
	//give the expected number of jumps (changes) from fromId to toId that occured along the specified branh length, 
	//in which the starting character is terminalStart and the terminal character is terminalEnd
	MDOUBLE getExpectation(const string& nodeName, int terminalStart, int terminalEnd, int fromId, int toId);
	//same as above, except here we return the probability of a jump from fromId to toId given 
	//terminal states terminalStart, terminalEnd in this branch
	MDOUBLE getProb(const string& nodeName, int terminalStart, int terminalEnd, int fromId, int toId);
    	
private:
	void init();
	void runOneIter(int state);
	void computeExpectationsAndPosterior();
	

private:

	//_node2Jumps: maps a node name (which specify a branch length) to 
	//the expected number of synonymous and nonsynonymous jumps between any two characters along the branch leading from the father to this node
	//given the terminal characters of this branch.
	//We use a "combined alphabet" to make access easier. see getCombinedState() for details
	//The dimension of the vector is the combined terminal state and the pair elements are: synonymous and non-synonymous jumps, respectively.

	map<string, vector<pair<MDOUBLE,MDOUBLE> > > _nodes2JumpsExp;
	
	//_node2JumpsProb: maps a node name (which specify a branch length) to 
	//the probability of a synonymous and non-synonymous jump between any two characters along the branch leading from the father to this node
	//given the terminal characters of this branch.
	//We use a "combined alphabet" to make access easier. see getCombinedState() for details
	//The dimension of the vector is the combined terminal state and the pair elements are: synonymous and non-synonymous jumps, respectively
	map<string, vector<pair<MDOUBLE,MDOUBLE> > > _nodes2JumpsProb;

};

#endif
