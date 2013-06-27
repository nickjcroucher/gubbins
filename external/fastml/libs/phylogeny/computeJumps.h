#ifndef ___COMPUTE_JUMPS__
#define ___COMPUTE_JUMPS__

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "alphabet.h"

#include <map>
#include <vector>
using namespace std;

/******************************************************************
This class compute jumps (events) by Suchard equations along differing branch lengths (according to a 
given tree), with the aim of giving the expectation of the number of jumps
from state a to state b given that the terminal states at the end of the branch are
x and y.
*******************************************************************/

class computeJumps  {
public:
	computeJumps(const MDOUBLE Lambda1, const MDOUBLE Lambda2, const MDOUBLE r=1);
	virtual ~computeJumps();

	//////////////////////////////////////////////////////////////////////////
	class gFunc  {
	public:
		gFunc(const MDOUBLE Lambda1, const MDOUBLE Lambda2 , const MDOUBLE r);
		gFunc(){};
		~gFunc(){};

		MDOUBLE gFunc_dr(MDOUBLE BranchLength);
		MDOUBLE g1Func_dr(MDOUBLE BranchLength);
		MDOUBLE g2Func_dr(MDOUBLE BranchLength);

		MDOUBLE g1Exp(MDOUBLE BranchLength);
		MDOUBLE g2Exp(MDOUBLE BranchLength);

	private:
		MDOUBLE _r;
		MDOUBLE _Lambda1;
		MDOUBLE _Lambda2;

		MDOUBLE _Alpha1;
		MDOUBLE _Alpha2;
		MDOUBLE _Alpha1_dr;
		MDOUBLE _Alpha2_dr;

		MDOUBLE _Alpha1_2;
		MDOUBLE _Alpha1_2_dr;

		MDOUBLE _delta;
		MDOUBLE _delta_dr;

		MDOUBLE _g1Part;
		MDOUBLE _g2Part;
		MDOUBLE _g1Part_dr;
		MDOUBLE _g2Part_dr;

	};
	//////////////////////////////////////////////////////////////////////////
	
	MDOUBLE getExpectation(const MDOUBLE BranchLength, int terminalStart, int terminalEnd, int fromId, int toId);
	MDOUBLE gainExp(MDOUBLE BranchLength,MDOUBLE prob01,MDOUBLE prob11);
	
	MDOUBLE gainExpGiven01(MDOUBLE BranchLength);
	MDOUBLE gainExpGiven00(MDOUBLE BranchLength);
	MDOUBLE gainExpGiven11(MDOUBLE BranchLength);
	MDOUBLE gainExpGiven10(MDOUBLE BranchLength);

	MDOUBLE lossExpGiven01(MDOUBLE BranchLength);
	MDOUBLE lossExpGiven00(MDOUBLE BranchLength);
	MDOUBLE lossExpGiven11(MDOUBLE BranchLength);
	MDOUBLE lossExpGiven10(MDOUBLE BranchLength);


	MDOUBLE gFunc_dr(MDOUBLE BranchLength);
    	
private:
	MDOUBLE m01(MDOUBLE BranchLength);
	MDOUBLE m00(MDOUBLE BranchLength);
	MDOUBLE m11(MDOUBLE BranchLength);
	MDOUBLE m10(MDOUBLE BranchLength);

	
	//MDOUBLE g1Func_dr(MDOUBLE BranchLength);
	//MDOUBLE g2Func_dr(MDOUBLE BranchLength);
	//MDOUBLE g1Exp(MDOUBLE BranchLength);
	//MDOUBLE g2Exp(MDOUBLE BranchLength);

	MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d);


	MDOUBLE _Lambda1;
	MDOUBLE _Lambda2;
	gFunc _gFuncStart0;
	gFunc _gFuncStart0MinusR;
	gFunc _gFuncStart1;
	gFunc _gFuncStart1MinusR;

};

#endif
