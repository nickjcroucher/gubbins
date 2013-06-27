#include "computePosteriorExpectationOfSubstitutions.h"
#include "definitions.h"
#include "computeDownAlg.h"
#include "computeUpAlg.h"
#include "matrixUtils.h"
#include "treeIt.h"
#include "likelihoodComputation.h"

using namespace std;

/********************************************************************************************
computePosteriorExpectationOfSubstitutions
*********************************************************************************************/
computePosteriorExpectationOfSubstitutions::computePosteriorExpectationOfSubstitutions(const tree &tr, const sequenceContainer &sc, const stochasticProcess *sp):
_tr(tr), _sc(sc){
	if(!sp){
		errorMsg::reportError("error in the constructor computePosteriorExpectationOfSubstitutions sp argument is NULL");
	}
	else{
		_sp = sp;
	}
}
/********************************************************************************************
Expectation of number of substitutions from character u to v --- =
sum over all substitutions x,y:
Posterior(Node=x,Father=y|D)*Exp(substitutions u to v|Node=x,Father=y)
The second term is given to the function as input (can be obtained via simulations)
*********************************************************************************************/
VVdouble computePosteriorExpectationOfSubstitutions::computeExpectationAcrossTree(
	simulateJumpsAbstract &sim,  //input given from simulation studies
	const VVVdouble &posteriorProbs,
	VVVdouble &expForBranch)
{
	//int numNodes = _tr.getNodesNum();
	int alphabetSize = _sp->alphabetSize();
	VVdouble res;
	resizeMatrix(res,alphabetSize,alphabetSize);
	treeIterTopDownConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		for (int fromState=0;fromState<alphabetSize;++fromState)
		{
			for (int toState=0;toState<alphabetSize;++toState)
			{
				if (fromState==toState) 
					continue;
				expForBranch[mynode->id()][fromState][toState] = computeExpectationOfChangePerBranch(sim,posteriorProbs,mynode,fromState,toState);
				res[fromState][toState] +=expForBranch[mynode->id()][fromState][toState];

			}
		}		
	}
	return res;
}
/********************************************************************************************
Posterior probabilities computed across entire tree, for all substitutions from character u to v 
*********************************************************************************************/
VVdouble computePosteriorExpectationOfSubstitutions::computePosteriorAcrossTree(
	simulateJumpsAbstract &sim, //input given from simulation studies
	const VVVdouble &posteriorProbsGivenTerminals,VVVdouble &probsForBranch)
{
	//int numNodes = _tr.getNodesNum();
	int alphabetSize = _sp->alphabetSize();
	// N: resized before 
	//probsForBranch.resize(numNodes);
	//for (int n=0;n<numNodes;++n)
	//	resizeMatrix(probsForBranch[n],alphabetSize,alphabetSize);

	VVdouble res;
	resizeMatrix(res,alphabetSize,alphabetSize);
	treeIterTopDownConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		for (int fromState=0;fromState<alphabetSize;++fromState)
		{
			for (int toState=0;toState<alphabetSize;++toState)
			{
				if (fromState==toState) 
					continue;
				probsForBranch[mynode->id()][fromState][toState]= computePosteriorOfChangePerBranch(sim,posteriorProbsGivenTerminals,mynode,fromState,toState);
				res[fromState][toState] +=probsForBranch[mynode->id()][fromState][toState];

			}
		}
	}
	return res;
}
/********************************************************************************************
*********************************************************************************************/
MDOUBLE computePosteriorExpectationOfSubstitutions::computePosteriorOfChangePerBranch(simulateJumpsAbstract &sim, //input given from simulation studies
	const VVVdouble &posteriorProbs,
	tree::nodeP node,
	int fromState, int toState)
{
	int alphabetSize = _sp->alphabetSize();
	MDOUBLE res = 0;

	for (int x=0;x<alphabetSize;++x)
	{
		for (int y=0;y<alphabetSize;++y)
		{
			res+=sim.getProb(node->name(),x,y,fromState,toState)*posteriorProbs[node->id()][x][y];
		}
	}
	return res;
}

/********************************************************************************************
Posterior of observing a certain state substitution along a branch:
P(Node=x,Father=y|D) = P(D,Node=x,Father=y)/P(D)
usage: posteriorPerNodePer2States[mynode->id()][fatherState][sonState]
*********************************************************************************************/
void computePosteriorExpectationOfSubstitutions::computePosteriorOfChangeGivenTerminals(VVVdouble &posteriorPerNodePer2States, int pos){
	int numNodes = _tr.getNodesNum();
	int alphabetSize = _sp->alphabetSize();
	posteriorPerNodePer2States.resize(numNodes);
	for (int n=0;n<posteriorPerNodePer2States.size();++n)
		resizeMatrix(posteriorPerNodePer2States[n],alphabetSize,alphabetSize);
	suffStatGlobalHomPos sscUp;
	suffStatGlobalHomPos sscDown; //for a reversible model
	sscUp.allocatePlace(numNodes,alphabetSize);
	computePijHom pi;
	pi.fillPij(_tr,*_sp); 

	computeUpAlg comp_Up;
	computeDownAlg comp_Down;
	comp_Up.fillComputeUp(_tr,_sc,pos,pi,sscUp);
	comp_Down.fillComputeDown(_tr,_sc,pos,pi,sscDown,sscUp);
	treeIterTopDownConst tIt(_tr);
	MDOUBLE ll = convert(likelihoodComputation::getLofPos(pos,_tr,_sc,pi,*_sp));
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		for (int sonState = 0; sonState<alphabetSize; ++sonState){
			for (int fatherState = 0; fatherState<alphabetSize; ++fatherState){
				posteriorPerNodePer2States[mynode->id()][fatherState][sonState]= computePosterioGivenTerminalsPerBranch(mynode->id(),sonState,fatherState,sscUp,sscDown, pi,ll,mynode->name());
			}
		}
	}
}
/********************************************************************************************
Posterior of observing a certain state substitution along a branch:
P(Node=sonState,Father=fatherState|D) = P(D,Node=sonState,Father=fatherState)/P(D)
usage: posteriorPerNodePer2States[mynode->id()][fatherState][sonState]
*********************************************************************************************/
MDOUBLE computePosteriorExpectationOfSubstitutions::computePosterioGivenTerminalsPerBranch
	(int nodeId,int sonState, int fatherState,suffStatGlobalHomPos &sscUp,
	suffStatGlobalHomPos &sscDown,computePijHom &pi, MDOUBLE &LLData, const string nodeName)
{
	MDOUBLE res, Down, Up, pij;
	Down = convert(sscDown.get(nodeId,fatherState));
	Up = convert(sscUp.get(nodeId,sonState));
	pij = pi.getPij(nodeId,fatherState,sonState);
	res=_sp->freq(fatherState)*Down*Up*pij;
	res/=LLData;
//	if(gainLossOptions::_printDEBUGinfo)
//		LOG(3,<<nodeName<<" son "<<sonState<<" Down "<<Down<<" father "<<fatherState<<" Up "<<Up<<" pij "<<pij<<" resDXY "<<resDXY<<" LLData "<<LLData<<" prob "<<res<<endl);

	if (res > 1 + 1e-4){
		LOGnOUT(3,<<nodeId<<" son "<<sonState<<" Down "<<Down<<" father "<<fatherState<<" Up "<<Up<<" pij "<<pij<<" res "<<res<<" LLData "<<LLData<<endl);
		res = 1;
	}
	if (res<-1e-4){
		LOGnOUT(3,<<nodeId<<" son "<<sonState<<" Down "<<Down<<" father "<<fatherState<<" Up "<<Up<<" pij "<<pij<<" res "<<res<<" LLData "<<LLData<<endl);
		res = 0;
	}
	if ((res > 1 + 0.000001) || (res<-0.000001)){
		string err = "Error in computePosteriorExpectationOfSubstitutions::computePosterioGivenTerminalsPerBranch, non probability value ";
		err+=double2string(res);
		err+=" at node ";
		err+=int2string(nodeId);
		err+=  " sonState ";
		err+= int2string(sonState);
		err+= " fatherState ";
		err+= int2string(fatherState);
		errorMsg::reportError(err);
	}
	return res;
}
/********************************************************************************************
*********************************************************************************************/
MDOUBLE computePosteriorExpectationOfSubstitutions::computeExpectationOfChangePerBranch(
	simulateJumpsAbstract &sim, //input given from simulation studies
	const VVVdouble &posteriorProbsGivenTerminals,
	tree::nodeP node,int fromState, int toState)
{
	int alphabetSize = _sp->alphabetSize();


	MDOUBLE nodeExpectation = 0;
	for (int x = 0; x<alphabetSize; ++x){
		for (int y = 0; y<alphabetSize; ++y){
			nodeExpectation+=(posteriorProbsGivenTerminals[node->id()][x][y]*
				sim.getExpectation(node->name(),x,y,fromState,toState));
			//DEBUG    
			LOG(6,<<"node "<<node->id()<<endl);
			LOG(6,<<"from "<<fromState<<" to "<<toState<<" given "<<x<<" and "<<y
				<<" post= "<<posteriorProbsGivenTerminals[node->id()][x][y]<<" sim= "<< sim.getExpectation(node->name(),x,y,fromState,toState)<<endl);
		}
	}
	return nodeExpectation;
}




