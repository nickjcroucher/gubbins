#ifndef ___SIMULATE_JUMPS_ABSTRACT_
#define ___SIMULATE_JUMPS_ABSTRACT_

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "alphabet.h"

#include <map>
#include <vector>
using namespace std;

/******************************************************************
This is an abstract class to various implementations of simulateJumps.
It was created to be a father class to the generic (original) implementation of
simulateJumps class simulateJumps (working on alphabets of either 0,1,2 or 0,1 
and class simulateCodonsJumps which is a variant simulateJumps that can handle the
61 sized alphabet without memory limitations.

The simulateJumps algorithm simulates jumps (events) along differing branch lengths (according to a 
given tree), with the aim of giving the expectation of the number of jumps
from state a to state b given that the terminal states at the end of the branch are
x and y.
*******************************************************************/

class simulateJumpsAbstract  {
public:
	simulateJumpsAbstract(const tree& inTree, const stochasticProcess& sp, const int alphabetSize);
	virtual ~simulateJumpsAbstract(){}
	virtual void runSimulation(int iterNum = 10000); 
	
	//for a branch length specified by a nodeName: 
	//give the expected number of jumps (changes) from fromId to toId that occured along the specified branh length, 
	//in which the starting character is terminalStart and the terminal character is terminalEnd
	virtual MDOUBLE getExpectation(const string& nodeName, int terminalStart, int terminalEnd, int fromId, int toId) = 0;
	//same as above, except here we return the probability of a jump from fromId to toId given 
	//terminal states terminalStart, terminalEnd in this branch
	virtual MDOUBLE getProb(const string& nodeName, int terminalStart, int terminalEnd, int fromId, int toId) = 0;
    	
protected:
	virtual int getCombinedState(int terminalStart, int terminalEnd) const;
	virtual int getCombinedAlphabetSize() const {return _alphabetSize*_alphabetSize;}
	virtual int getStartId(int combinedState) const;
	virtual int getEndId(int combinedState) const;

	virtual void init() = 0;
	virtual void runOneIter(int state) = 0;
	virtual void computeExpectationsAndPosterior() = 0;

	// a comparison function to be used in sort init
	static bool compareDist(tree::nodeP node1, tree::nodeP node2){	return (node1->dis2father() < node2->dis2father());}
	

protected:
	tree _tree;
	stochasticProcess _sp;
	const int _alphabetSize;

	Vdouble _waitingTimeParams;//each entry is the lambda parameter of the exponential distribution modeling the waiting time for "getting out" of state i

	//_jumpProbs[i][j] is the probability of jumping from state i to state j (given that a change has ocured).
	VVdouble _jumpProbs; 

	//the number of times we reached a certain combination of terminal states for each branch lengths
	//e.g. the number of times we observed 0,1 at terminal states given branch length 0.03
	//this is used to to afterwards normalize (i.e. compute the expectation) the _nodes2JumpsExp values
	map<string, Vdouble> _totalTerminals; 

	vector<tree::nodeP> _orderNodesVec; //internal use: the branch are sorted in ascending order 

};

#endif
