// $Id: bblEM.cpp 6110 2009-04-26 14:06:02Z cohenofi $
#include "bblEM.h"
#include "likelihoodComputation.h"
using namespace likelihoodComputation;
#include "computeUpAlg.h"
#include "computeDownAlg.h"
#include "computeCounts.h"
#include "treeIt.h"
#include "fromCountTableComponentToDistance.h"
#include <ctime>

bblEM::bblEM(tree& et,
			 const sequenceContainer& sc,
			 const stochasticProcess& sp,
			 const Vdouble * weights,
			 const int maxIterations,
			 const MDOUBLE epsilon,
			 const MDOUBLE tollForPairwiseDist,
			 unObservableData*  _unObservableData_p) :
_et(et),_sc(sc),_sp(sp),_weights(weights),_unObservableData_p(_unObservableData_p) 
{
	_treeLikelihood = compute_bblEM(maxIterations,epsilon,tollForPairwiseDist);
}


MDOUBLE bblEM::compute_bblEM(
			const int maxIterations,
			const MDOUBLE epsilon,
			const MDOUBLE tollForPairwiseDist){
	allocatePlace();
	MDOUBLE oldL=VERYSMALL;
	MDOUBLE currL = VERYSMALL;
	tree oldT = _et;
	for (int i=0; i < maxIterations; ++i) {
		computeUp();
		currL = likelihoodComputation::getTreeLikelihoodFromUp2(_et,_sc,_sp,_cup,_posLike,_weights,_unObservableData_p);
		//////////////
		//MDOUBLE checkUpLL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et, _sc, _sp, _weights);
		//cerr << "checkUpLL = "<<checkUpLL <<" curll = "<<currL<<endl;
		///////////////
		oldT = _et;
		if (currL < oldL + epsilon) { // need to break
			if (currL<oldL) {
				_et = oldT;
				if(_unObservableData_p)
					_unObservableData_p->setLforMissingData(_et,&_sp);
				return oldL; // keep the old tree, and old likelihood
			} else {
                //update the tree and likelihood and return
				return currL;
			}
		}
		bblEM_it(tollForPairwiseDist);
		oldL = currL;
	}
	// in the case were we reached max_iter, we have to recompute the likelihood of the new tree...
	computeUp();
	if(_unObservableData_p){
        _unObservableData_p->setLforMissingData(_et,&_sp);
	}
	currL = likelihoodComputation::getTreeLikelihoodFromUp2(_et,_sc,_sp,_cup,_posLike,_weights,_unObservableData_p);
	if (currL<oldL) {
		_et = oldT;
		if(_unObservableData_p)
			_unObservableData_p->setLforMissingData(_et,&_sp);
		return oldL; // keep the old tree, and old likelihood
	} 
	else 
        return currL;
}

void bblEM::allocatePlace() {
	_computeCountsV.resize(_et.getNodesNum()); //initiateTablesOfCounts
	for (int i=0; i < _computeCountsV.size(); ++i) {
		_computeCountsV[i].countTableComponentAllocatePlace(_sp.alphabetSize(),_sp.categories());
	}
	_cup.allocatePlace(_sc.seqLen(),_sp.categories(), _et.getNodesNum(), _sc.alphabetSize());
	_cdown.allocatePlace(_sp.categories(), _et.getNodesNum(), _sc.alphabetSize());
}

void bblEM::bblEM_it(const MDOUBLE tollForPairwiseDist){
	for (int i=0; i < _computeCountsV.size(); ++i) {
		_computeCountsV[i].zero();
	}
	for (int i=0; i < _sc.seqLen(); ++i) {
		computeDown(i);
		addCounts(i); // computes the counts and adds to the table.
	}
	optimizeBranches(tollForPairwiseDist);
	if(_unObservableData_p){
		_unObservableData_p->setLforMissingData(_et,&_sp);
	}
}

void bblEM::optimizeBranches(const MDOUBLE tollForPairwiseDist){
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			fromCountTableComponentToDistance from1(_computeCountsV[mynode->id()],_sp,tollForPairwiseDist,mynode->dis2father(),_unObservableData_p);
			from1.computeDistance();
			mynode->setDisToFather(from1.getDistance());
			if(_unObservableData_p){
				_unObservableData_p->setLforMissingData(_et,&_sp);
			}
		}
	}
}

void bblEM::computeUp(){
	_pij.fillPij(_et,_sp,0); // 0 is becaues we compute Pij(t) and not its derivations...
	computeUpAlg cupAlg;
	for (int pos=0; pos < _sc.seqLen(); ++pos) {
        for (int categor = 0; categor < _sp.categories(); ++categor) {
			cupAlg.fillComputeUp(_et,_sc,pos,_pij[categor],_cup[pos][categor]);
		}
	}
 }

void bblEM::computeDown(const int pos){
	computeDownAlg cdownAlg;
	for (int categor = 0; categor < _sp.categories(); ++categor) {
		cdownAlg.fillComputeDown(_et,_sc,pos,_pij[categor],_cdown[categor],_cup[pos][categor]);
	}
 }

void bblEM::addCounts(const int pos){
	//MDOUBLE posProb = 
	//	likelihoodComputation::getProbOfPosWhenUpIsFilledGam(pos,_et,_sc,_sp,_cup);
						
	MDOUBLE weig = (_weights ? (*_weights)[pos] : 1.0);
	if (weig == 0) return;
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			addCounts(pos,mynode,_posLike[pos],weig);
		}
	}
}

void bblEM::addCounts(const int pos, tree::nodeP mynode, const doubleRep posProb, const MDOUBLE weig){

	computeCounts cc;
	for (int categor =0; categor< _sp.categories(); ++ categor) {
			cc.computeCountsNodeFatherNodeSonHomPos(_sc,
										_pij[categor],
										_sp,
										_cup[pos][categor],
										_cdown[categor],
										weig,
										posProb,
										mynode,
										_computeCountsV[mynode->id()][categor],
										_sp.ratesProb(categor));
	}
}          

