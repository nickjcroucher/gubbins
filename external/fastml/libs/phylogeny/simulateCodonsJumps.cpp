#include "simulateCodonsJumps.h"
#include "talRandom.h"
#include "someUtil.h"
#include "codon.h"
#include <algorithm>


simulateCodonsJumps::simulateCodonsJumps(const tree& inTree, const stochasticProcess& sp, const int alphabetSize)
: simulateJumpsAbstract(inTree,sp,alphabetSize)	
{
}

simulateCodonsJumps::~simulateCodonsJumps()
{
}

void simulateCodonsJumps::init()
{
	//init the vector of waiting times. 
	_waitingTimeParams.clear();
	_waitingTimeParams.resize(_alphabetSize);
	int i, j;
	for (i = 0; i < _alphabetSize; ++i)
	{
		_waitingTimeParams[i] = -_sp.dPij_dt(i, i, 0.0);
		
	}

	//init _jumpProbs.
	_jumpProbs.clear();
	_jumpProbs.resize(_alphabetSize);
	for (i = 0; i < _alphabetSize; ++i)
	{
		MDOUBLE sum = 0.0;
		_jumpProbs[i].resize(_alphabetSize);
		for (j = 0; j < _alphabetSize; ++j)
		{
			if (i == j)
				_jumpProbs[i][j] = 0.0;
			else
			{
				_jumpProbs[i][j] = _sp.dPij_dt(i, j, 0.0) / _waitingTimeParams[i];
			}
			sum += _jumpProbs[i][j];
		}
		if (! DEQUAL(sum, 1.0)){
			string err = "error in simulateJumps::init(): sum probabilities is not 1 and equal to ";
			err+=double2string(sum);
			errorMsg::reportError(err);
		}
	}

	//init _orderNodesVec: a vector in which the branch lengths are ordered in ascending order
	_tree.getAllNodes(_orderNodesVec, _tree.getRoot());
	sort(_orderNodesVec.begin(), _orderNodesVec.end(), simulateJumpsAbstract::compareDist); 

	_nodes2JumpsExp.clear();
	_nodes2JumpsProb.clear();
//
	vector<pair<MDOUBLE,MDOUBLE> > zeroCombinedStates2jumps;
	for(i = 0;i < getCombinedAlphabetSize();++i){
		pair<MDOUBLE,MDOUBLE> syn_and_nonSyn_jumps(0.0,0.0);
		zeroCombinedStates2jumps.push_back(syn_and_nonSyn_jumps);
	}
	Vdouble zeroVector(getCombinedAlphabetSize(),0.0);
	for (i = 0; i < _orderNodesVec.size(); ++i)
	{
		string nodeName = _orderNodesVec[i]->name();
		_nodes2JumpsExp[nodeName] = zeroCombinedStates2jumps;
		_nodes2JumpsProb[nodeName] = zeroCombinedStates2jumps;
		for (j=0; j<getCombinedAlphabetSize();++j)
			_totalTerminals[nodeName]=zeroVector;
	}		
}


//simulate jumps starting from startState. The simulation continue until the maxTime is reached. In each step:
//1. Draw a new waiting time.
//2. Go over all branches shorter than nextJumpTime and update their jumpsNum between the states that were switched 
//	(these branches will not be affected by the current jump): 
//	however they might have been affected by the previous jump
//3. Draw a new state
void simulateCodonsJumps::runOneIter(int startState)
{
	int substitutionType;
	MDOUBLE maxTime = _orderNodesVec[_orderNodesVec.size()-1]->dis2father();
	MDOUBLE totalTimeTillJump = 0.0;
	int curState = startState;
	int smallestBranchNotUpdatedSofar = 0;
	vector<pair<int, int> > jumpsSoFar(0);
	while (totalTimeTillJump < maxTime)
	{
		MDOUBLE avgWaitingTime = 1 / _waitingTimeParams[curState];
		MDOUBLE nextJumpTime = totalTimeTillJump + talRandom::rand_exp(avgWaitingTime);
		//go over all branches that "finished" their simulation (shorter than nextJumpTime) and update with their _nodes2JumpsExp 
		//with the jumps that occured between the terminal Ids: startState-->curState
		for (int b = smallestBranchNotUpdatedSofar; b < _orderNodesVec.size(); ++b)
		{
			if (_orderNodesVec[b]->dis2father() > nextJumpTime)
			{
				smallestBranchNotUpdatedSofar = b;
				break;
			}
			string nodeName = _orderNodesVec[b]->name();
			//update all the jumps that occured along the branch
			int terminalState = getCombinedState(startState, curState);
			_totalTerminals[nodeName][terminalState]++;
			//update all longer branches with all jumps that occurred till now
			vector<bool> jumpsSoFarBool(getCombinedAlphabetSize(),false);
			for (int j = 0; j < jumpsSoFar.size(); ++j)
			{
				int combinedJumpState = getCombinedState(jumpsSoFar[j].first, jumpsSoFar[j].second);
				jumpsSoFarBool[combinedJumpState]=true;
				if(substitutionType == 1)
					_nodes2JumpsExp[nodeName][terminalState].first += 1;	
				else if(substitutionType == 2)
					_nodes2JumpsExp[nodeName][terminalState].second += 1;
			}
			for (int combined=0;combined<jumpsSoFarBool.size();++combined)
			{
				if (jumpsSoFarBool[combined]){
					if(substitutionType == 1)
						_nodes2JumpsProb[nodeName][terminalState].first += 1;	
					else if(substitutionType == 2)
						_nodes2JumpsProb[nodeName][terminalState].second += 1;
				}
			}
		}
		totalTimeTillJump = nextJumpTime;
		int nextState = giveRandomState(_alphabetSize,curState,_jumpProbs);
		substitutionType = codonUtility::codonReplacement(curState,nextState);
		jumpsSoFar.push_back(pair<int,int>(curState, nextState));
		curState = nextState;
	}
}


void simulateCodonsJumps::computeExpectationsAndPosterior(){
	//scale _nodes2JumpsExp so it will represent expectations
	map<string, vector<pair<MDOUBLE,MDOUBLE> > >::iterator iterExp = _nodes2JumpsExp.begin();
	for (; iterExp != _nodes2JumpsExp.end(); ++iterExp)
	{//each node
		string nodeName = iterExp->first;
		for (int termState = 0; termState < getCombinedAlphabetSize(); ++termState)
		{
			map<string, Vdouble>::iterator iterTerm = _totalTerminals.find(nodeName);
			map<string, vector<pair<MDOUBLE,MDOUBLE> > >::iterator iterProb = _nodes2JumpsProb.find(nodeName);
			if ((iterTerm==_totalTerminals.end()) || (iterProb==_nodes2JumpsProb.end()))
			{
				errorMsg::reportError("error in simulateJumps::runSimulation, unknown reason: cannot find nodeName in map");
			}

			if (iterTerm->second[termState]==0){ //never reached these terminal states
				if((iterExp->second[termState].first == 0)&&(iterExp->second[termState].second == 0)&&
					((iterProb->second[termState].first == 0)&&(iterProb->second[termState].second == 0))) continue;
				else
					errorMsg::reportError("error in simulateCodonJumps::runSimulation, 0 times reached termState but non-zero for jumpCount");
			}
			(iterExp->second[termState].first) /= iterTerm->second[termState];
			(iterExp->second[termState].second) /= iterTerm->second[termState];
			(iterProb->second[termState].first) /= iterTerm->second[termState];
			(iterProb->second[termState].second) /= iterTerm->second[termState];
		}
	}
}


MDOUBLE simulateCodonsJumps::getExpectation(const string& nodeName, int terminalStart, int terminalEnd, int fromId, int toId)
{
	//map <string, VVdouble>::iterator pos;//Old
	map<string, vector<pair<MDOUBLE,MDOUBLE> > >::iterator pos;
	if ((pos = _nodes2JumpsExp.find(nodeName)) == _nodes2JumpsExp.end())
	{
		string err="error in simulateCodonJumps::getExpectation: cannot find node "+nodeName;
		errorMsg::reportError(err);
	}
	int combinedTerminalState = getCombinedState(terminalStart, terminalEnd);
	//Old
	//int combinedJumpState = getCombinedState(fromId, toId);
	//return (pos->second[combinedTerminalState][combinedJumpState]);

	MDOUBLE expectation;
	if(codonUtility::codonReplacement(fromId,toId) == 1)
		expectation = pos->second[combinedTerminalState].first;
	else if(codonUtility::codonReplacement(fromId,toId) == 2)
		expectation = pos->second[combinedTerminalState].second;
	return (expectation);
}


MDOUBLE simulateCodonsJumps::getProb(const string& nodeName, int terminalStart, int terminalEnd, int fromId, int toId){
	//map <string, VVdouble>::iterator pos;
	map<string, vector<pair<MDOUBLE,MDOUBLE> > >::iterator pos;
	if ((pos = _nodes2JumpsProb.find(nodeName)) == _nodes2JumpsProb.end())
	{
		string err="error in simulateCodonJumps::getProb: cannot find node "+nodeName;
		errorMsg::reportError(err);
	}
	int combinedTerminalState = getCombinedState(terminalStart, terminalEnd);
	//Old
	//int combinedJumpState = getCombinedState(fromId, toId);
	//return (pos->second[combinedTerminalState][combinedJumpState]);

	MDOUBLE expectation;
	if(codonUtility::codonReplacement(fromId,toId) == 1)
		expectation = pos->second[combinedTerminalState].first;
	else if(codonUtility::codonReplacement(fromId,toId) == 2)
		expectation = pos->second[combinedTerminalState].second;
	return (expectation);
}