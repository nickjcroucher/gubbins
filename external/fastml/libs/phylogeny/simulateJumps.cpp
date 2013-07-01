#include "simulateJumps.h"
#include "talRandom.h"
#include "someUtil.h"
#include <algorithm>


simulateJumps::simulateJumps(const tree& inTree, const stochasticProcess& sp, const int alphabetSize)
: simulateJumpsAbstract(inTree,sp,alphabetSize)	
{
}

simulateJumps::~simulateJumps()
{
}

void simulateJumps::init()
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
	//_jumpProbs[i][j] = Q[i][j] / -Q[i][i]
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
	VVdouble zeroMatrix(getCombinedAlphabetSize());
	for (i = 0; i < getCombinedAlphabetSize(); ++i)
		zeroMatrix[i].resize(getCombinedAlphabetSize(), 0.0);
	Vdouble zeroVector(getCombinedAlphabetSize(),0.0);
	for (i = 0; i < _orderNodesVec.size(); ++i)
	{
		string nodeName = _orderNodesVec[i]->name();
		_nodes2JumpsExp[nodeName] = zeroMatrix;
		_nodes2JumpsProb[nodeName] = zeroMatrix;
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
void simulateJumps::runOneIter(int startState)
{
	MDOUBLE maxTime = _orderNodesVec[_orderNodesVec.size()-1]->dis2father();
	MDOUBLE totalTimeTillJump = 0.0;
	int jumpsNum = 0;
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
                _nodes2JumpsExp[nodeName][terminalState][combinedJumpState] += 1;
			}
			for (int combined=0;combined<jumpsSoFarBool.size();++combined)
			{
				if (jumpsSoFarBool[combined])
					_nodes2JumpsProb[nodeName][terminalState][combined]+=1;
			}
		}
		totalTimeTillJump = nextJumpTime;
		int nextState = giveRandomState(_alphabetSize,curState, _jumpProbs);
		jumpsSoFar.push_back(pair<int,int>(curState, nextState));
		curState = nextState;
		++jumpsNum;
	}
}


void simulateJumps::computeExpectationsAndPosterior(){
	//scale _nodes2JumpsExp so it will represent expectations
	map<string, VVdouble>::iterator iterExp = _nodes2JumpsExp.begin();
	for (; iterExp != _nodes2JumpsExp.end(); ++iterExp)
	{
		string nodeName = iterExp->first;
		for (int termState = 0; termState < getCombinedAlphabetSize(); ++termState)
		{
			for (int jumpState = 0; jumpState < getCombinedAlphabetSize(); ++jumpState)
			{

				//(iter->second[termState][jumpState]) /= static_cast<MDOUBLE>(iterNum);
				map<string, Vdouble>::iterator iterTerm = _totalTerminals.find(nodeName);
				map<string, VVdouble>::iterator iterProb = _nodes2JumpsProb.find(nodeName);
				if ((iterTerm==_totalTerminals.end()) || (iterProb==_nodes2JumpsProb.end()))
				{
					errorMsg::reportError("error in simulateJumps::runSimulation, unknown reason: cannot find nodeName in map");
				}
				if ((iterTerm->second[termState]==0)){ //never reached these terminal states
					if ((iterExp->second[termState][jumpState]==0) && (iterProb->second[termState][jumpState]==0))
						continue;//leave the value of _nodes2JumpsExp and _nodes2JumpsProb as zero
					else {
						errorMsg::reportError("error in simulateJumps::runSimulation, 0 times reached termState but non-zero for jumpCount");
					}
				}
				(iterExp->second[termState][jumpState]) /= iterTerm->second[termState];
				
				(iterProb->second[termState][jumpState]) /= iterTerm->second[termState];
				
			}
		}
	}
}


MDOUBLE simulateJumps::getExpectation(const string& nodeName, int terminalStart, int terminalEnd, int fromId, int toId)
{
	map <string, VVdouble>::iterator pos;
	if ((pos = _nodes2JumpsExp.find(nodeName)) == _nodes2JumpsExp.end())
	{
		string err="error in simulateJumps::getExpectation: cannot find node "+nodeName;
		errorMsg::reportError(err);
	}
	int combinedTerminalState = getCombinedState(terminalStart, terminalEnd);
	int combinedJumpState = getCombinedState(fromId, toId);
	return (pos->second[combinedTerminalState][combinedJumpState]);
}


MDOUBLE simulateJumps::getProb(const string& nodeName, int terminalStart, int terminalEnd, int fromId, int toId){
	map <string, VVdouble>::iterator pos;
	if ((pos = _nodes2JumpsProb.find(nodeName)) == _nodes2JumpsProb.end())
	{
		string err="error in simulateJumps::getProb: cannot find node "+nodeName;
		errorMsg::reportError(err);
	}
	int combinedTerminalState = getCombinedState(terminalStart, terminalEnd);
	int combinedJumpState = getCombinedState(fromId, toId);
	return (pos->second[combinedTerminalState][combinedJumpState]);
}