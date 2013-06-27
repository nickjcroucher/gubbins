#include "definitions.h"
#include "computeDownAlg.h"
#include "computeUpAlg.h"
#include "matrixUtils.h"
#include "treeIt.h"
#include "likelihoodComputation.h"
#include "computePosteriorExpectationOfSubstitutions_nonReversibleSp.h"

using namespace std;



/********************************************************************************************
Posterior of observing a certain state substitution along a branch:
P(Node=x,Father=y|D) = P(D,Node=x,Father=y)/P(D)
usage: posteriorPerNodePer2States[mynode->id()][fatherState][sonState]
*********************************************************************************************/
void computePosteriorExpectationOfSubstitutions_nonReversibleSp::computePosteriorOfChangeGivenTerminals(VVVdouble &posteriorPerNodePer2States, int pos){
	int numNodes = _tr.getNodesNum();
	int alphabetSize = _sp->alphabetSize();
	posteriorPerNodePer2States.resize(numNodes);
	for (int n=0;n<posteriorPerNodePer2States.size();++n)
		resizeMatrix(posteriorPerNodePer2States[n],alphabetSize,alphabetSize);
	suffStatGlobalHomPos sscUp;
	suffStatGlobalGamPos sscDownNonRev;	// The "Gam" is used for the letter at father - sscGivenRoot
	sscUp.allocatePlace(numNodes,alphabetSize);
	computePijHom pi;
	pi.fillPij(_tr,*_sp); 

	computeUpAlg comp_Up;
	computeDownAlg comp_Down;
	comp_Up.fillComputeUp(_tr,_sc,pos,pi,sscUp);
	comp_Down.fillComputeDownNonReversible(_tr,_sc,pos,pi,sscDownNonRev,sscUp);
	treeIterTopDownConst tIt(_tr);
	MDOUBLE ll = convert(likelihoodComputation::getLofPos(pos,_tr,_sc,pi,*_sp));
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		for (int sonState = 0; sonState<alphabetSize; ++sonState){
			for (int fatherState = 0; fatherState<alphabetSize; ++fatherState){
				posteriorPerNodePer2States[mynode->id()][fatherState][sonState]= computePosterioGivenTerminalsPerBranch(mynode->id(),sonState,fatherState,sscUp,sscDownNonRev, pi,ll,mynode->name());
			}
		}
	}
}

/********************************************************************************************
Posterior of observing a certain state substitution along a branch:
P(Node=sonState,Father=fatherState|D) = P(D,Node=sonState,Father=fatherState)/P(D)
usage: posteriorPerNodePer2States[mynode->id()][fatherState][sonState]
*********************************************************************************************/
MDOUBLE computePosteriorExpectationOfSubstitutions_nonReversibleSp::computePosterioGivenTerminalsPerBranch
	(int nodeId,int sonState, int fatherState,suffStatGlobalHomPos &sscUp,
	suffStatGlobalGamPos &sscDown,computePijHom &pi, MDOUBLE &LLData, const string nodeName)
{
	MDOUBLE res=0.0;
	MDOUBLE resDXY, Down, Up, pij;
	for (int stateAtRoot = 0; stateAtRoot<_sp->alphabetSize(); ++stateAtRoot){
		Down = convert(sscDown.get(stateAtRoot,nodeId,fatherState));
		Up = convert(sscUp.get(nodeId,sonState));
		pij = pi.getPij(nodeId,fatherState,sonState);

		res+=(_sp->freq(stateAtRoot)*
			Down*
			Up*
			pij);
	}
	resDXY = res;
	res/=LLData;
//	if(gainLossOptions::_printDEBUGinfo)
//		LOG(3,<<nodeName<<" son "<<sonState<<" Down "<<Down<<" father "<<fatherState<<" Up "<<Up<<" pij "<<pij<<" resDXY "<<resDXY<<" LLData "<<LLData<<" prob "<<res<<endl);

	if (res > 1 + 1e-4){
		LOGnOUT(3,<<nodeId<<" son "<<sonState<<" Down "<<Down<<" father "<<fatherState<<" Up "<<Up<<" pij "<<pij<<" resDXY "<<resDXY<<" LLData "<<LLData<<" prob "<<res<<endl);
		res = 1;
	}
	if (res<-1e-4){
		LOGnOUT(3,<<nodeId<<" son "<<sonState<<" Down "<<Down<<" father "<<fatherState<<" Up "<<Up<<" pij "<<pij<<" resDXY "<<resDXY<<" LLData "<<LLData<<" prob "<<res<<endl);
		res = 0;
	}
	if ((res > 1 + 0.000001) || (res<-0.000001)){
		string err = "Error in computePosteriorExpectationOfSubstitutions_nonReversibleSp::computePosterioGivenTerminalsPerBranch, non probability value ";
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