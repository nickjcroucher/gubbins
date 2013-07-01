// $Id: Nni.cpp 962 2006-11-07 15:13:34Z privmane $

// version 1.00
// last modified 3 Nov 2002
#include "definitions.h"
#include "treeUtil.h"
#include "treeIt.h"
#include "Nni.h"
#include "bblEM.h"
#include "logFile.h"
#include <algorithm>
#include <iostream>
using namespace std;

NNI::NNI(const sequenceContainer& sc,
				   const stochasticProcess& sp,
				const Vdouble * weights): _sc(sc),_sp(sp),_weights(weights) {
	_bestScore = VERYSMALL;
}


tree NNI::NNIstep(tree et) {
	et.create_names_to_internal_nodes();
	treeIterTopDown tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (mynode->isLeaf() || mynode->isRoot()) continue; // swaping only internal nodes
		tree newT1 = NNIswap1(et,mynode);
		tree newT2 = NNIswap2(et,mynode);
		MDOUBLE treeScore1 = evalTree(newT1,_sc);
		MDOUBLE treeScore2 = evalTree(newT2,_sc);
		if (treeScore1 > _bestScore) {
			_bestTree = newT1;
			_bestScore = treeScore1;
			LOG(5,<<"new Best Tree: "<<_bestScore<<endl);
			LOGDO(5,et.output(myLog::LogFile()));
		}
		if (treeScore2 > _bestScore) {
			_bestTree = newT2;
			_bestScore = treeScore2;
			LOG(5,<<"new Best Tree: "<<_bestScore<<endl);
			LOGDO(5,et.output(myLog::LogFile()));
		}
	}
	return _bestTree;
}

tree NNI::NNIswap1(tree et,tree::nodeP mynode) {
	tree::nodeP mynodeInNewTree = et.findNodeByName(mynode->name());
#ifdef VERBOS
	LOG(5,<<"b4 swap1"<<endl);
	LOGDO(5,et.output(myLog::LogFile()));
#endif

	tree::nodeP fatherNode = mynodeInNewTree->father();
	tree::nodeP nodeToSwap1 = mynodeInNewTree->father()->getSon(0);
	// it might be me
	if (nodeToSwap1 == mynodeInNewTree) 
		nodeToSwap1 = mynodeInNewTree->father()->getSon(1);
	tree::nodeP nodeToSwap2 = mynodeInNewTree->getSon(0);

	et.removeNodeFromSonListOfItsFather(nodeToSwap1);
	et.removeNodeFromSonListOfItsFather(nodeToSwap2);
	nodeToSwap2->setFather(fatherNode);
	fatherNode->setSon(nodeToSwap2);
	nodeToSwap1->setFather(mynodeInNewTree);
	mynodeInNewTree->setSon(nodeToSwap1);
#ifdef VERBOS
	LOG(5,<<"after swap1"<<endl);
	LOGDO(5,et.output(myLog::LogFile()));
#endif
	
	return et;
}

tree NNI::NNIswap2(tree et,tree::nodeP mynode) {
#ifdef VERBOS
	LOG(5,<<"b4 swap2"<<endl);
	LOGDO(5,et.output(myLog::LogFile()));
#endif
	tree::nodeP mynodeInNewTree = et.findNodeByName(mynode->name());


	tree::nodeP fatherNode = mynodeInNewTree->father();
	tree::nodeP nodeToSwap1 = mynodeInNewTree->father()->getSon(0);
	// it might be me
	if (nodeToSwap1 == mynodeInNewTree) 
		nodeToSwap1 = mynodeInNewTree->father()->getSon(1);
	tree::nodeP nodeToSwap2 = mynodeInNewTree->getSon(1);
	et.removeNodeFromSonListOfItsFather(nodeToSwap1);
	et.removeNodeFromSonListOfItsFather(nodeToSwap2);
	nodeToSwap2->setFather(fatherNode);
	fatherNode->setSon(nodeToSwap2);
	nodeToSwap1->setFather(mynodeInNewTree);
	mynodeInNewTree->setSon(nodeToSwap1);
#ifdef VERBOS
	LOG(5,<<"after swap2"<<endl);
	LOGDO(5,et.output(myLog::LogFile()));
#endif //VERBOS
	return et;

}





MDOUBLE NNI::evalTree(tree& et,const sequenceContainer& sc) {
#ifdef VERBOS
	LOG(5,<<"b4 bbl in alltrees"<<endl);
	LOGDO(5,et.output(myLog::LogFile()));
#endif
	bblEM bblEM1(et,sc,_sp,_weights);
	MDOUBLE res = bblEM1.getTreeLikelihood();
	return res;
}	




