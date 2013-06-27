// $Id: NNiProp.cpp 962 2006-11-07 15:13:34Z privmane $

#include "definitions.h"
#include "treeIt.h"
#include "treeUtil.h"
#include "NNiProp.h"
#include "bblEM.h"
#include "bblEMProportional.h"
#include "logFile.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
using namespace std;

NNiProp::NNiProp(vector<sequenceContainer>& sc,
					 vector<stochasticProcess>& sp,
					const vector<Vdouble *> * weights,
					vector<char>* nodeNotToSwap):_nodeNotToSwap(nodeNotToSwap),
						_sc(sc),_sp(sp),_weights(weights) {
	_bestScore = VERYSMALL;
	_treeEvaluated =-1;
	_out = NULL;
	
}

void NNiProp::setOfstream(ostream* out) {
	_out = out;
}

tree NNiProp::NNIstep(tree et) {
	et.create_names_to_internal_nodes();
	_bestScore = evalTree(et);
	_bestTree = et;
	treeIterTopDown tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (mynode->isLeaf() || mynode->isRoot()) continue; // swaping only internal nodes

		if (_nodeNotToSwap) {
			if ((*_nodeNotToSwap)[mynode->id()]) {
				continue;
			}
		}
		tree newT1 = NNIswap1(et,mynode);
		tree newT2 = NNIswap2(et,mynode);
		MDOUBLE treeScore1 = evalTree(newT1);
		if (treeScore1 > _bestScore) {
			_bestTree = newT1;
			_bestScore = treeScore1;
			LOG(5,<<"new Best Tree: "<<_bestScore<<endl);
			if (_out) (*_out)<<"new Best Tree: "<<_bestScore<<endl;
			_bestTree.output(*_out);

		}
		MDOUBLE treeScore2 = evalTree(newT2);
		if (treeScore2 > _bestScore) {
			_bestTree = newT2;
			_bestScore = treeScore2;
			LOG(5,<<"new Best Tree: "<<_bestScore<<endl);
			if (_out) (*_out)<<"new Best Tree: "<<_bestScore<<endl;
			_bestTree.output(*_out);
		}
	}
	return _bestTree;
}

tree NNiProp::NNIswap1(tree et,tree::nodeP mynode) {
	tree::nodeP mynodeInNewTree = et.findNodeByName(mynode->name());
#ifdef VERBOS
	LOG(5,<<"b4 swap1"<<endl);
	LOGDO(5,et.output(myLog::LogFile()));
#endif

	tree::nodeP fatherNode = mynodeInNewTree->father();
	tree::nodeP nodeToSwap1 = mynodeInNewTree->father()->getSon(0);
	// it might be me
	if (nodeToSwap1 == mynodeInNewTree) nodeToSwap1 = mynodeInNewTree->father()->getSon(1);
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

tree NNiProp::NNIswap2(tree et,tree::nodeP mynode) {
#ifdef VERBOS
	LOG(5,<<"b4 swap2"<<endl);
	LOGDO(5,et.output(myLog::LogFile()));
#endif
	tree::nodeP mynodeInNewTree = et.findNodeByName(mynode->name());


	tree::nodeP fatherNode = mynodeInNewTree->father();
	tree::nodeP nodeToSwap1 = mynodeInNewTree->father()->getSon(0);
	// it might be me
	if (nodeToSwap1 == mynodeInNewTree) nodeToSwap1 = mynodeInNewTree->father()->getSon(1);
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
#endif
	return et;

}

MDOUBLE NNiProp::evalTree(tree& et) {
#ifdef VERBOS
	LOG(5,<<"b4 bbl in alltrees"<<endl);
	LOGDO(5,et.output(myLog::LogFile()));
#endif
	bblEMProportional bblEMprop1(et,_sc,_sp,_weights);
	MDOUBLE res = bblEMprop1.getTreeLikelihood();
//	MDOUBLE res = 12;
	_treeEvaluated++;
	//	cerr.precision(5);
	_out->precision(5);

	if (_treeEvaluated) LOG(5,<<"tree: "<<_treeEvaluated<< "score = "<<res<<endl);
	if ((_out)&&(_treeEvaluated)) (*_out)<<"tree: "<<_treeEvaluated<< "score = "<<res<<endl;
	return res;
}	




