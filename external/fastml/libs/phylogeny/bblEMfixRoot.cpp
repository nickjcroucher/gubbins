// $Id: bblEM.cpp 4478 2008-07-17 17:09:55Z cohenofi $
#include "bblEMfixRoot.h"
#include "likelihoodComputation.h"
using namespace likelihoodComputation;
#include "computeUpAlg.h"
#include "computeDownAlg.h"
#include "computeCounts.h"
#include "treeIt.h"
#include "fromCountTableComponentToDistancefixRoot.h"
#include <ctime>

bblEMfixRoot::bblEMfixRoot(tree& et,
				const sequenceContainer& sc,
				const stochasticProcess& sp,
				const Vdouble * weights,
				const int maxIterations,
				const MDOUBLE epsilon,
				const MDOUBLE tollForPairwiseDist,
				unObservableData*  unObservableData_p) :
_et(et),_sc(sc),_sp(sp),_weights (weights),_unObservableData_p(unObservableData_p) 
{
	//if(!plogLforMissingData){
	//	_plogLforMissingData = NULL;
	//}
	_treeLikelihood = compute_bblEM(maxIterations,epsilon,tollForPairwiseDist);
}


MDOUBLE bblEMfixRoot::compute_bblEM(
			const int maxIterations,
			const MDOUBLE epsilon,
			const MDOUBLE tollForPairwiseDist){
	allocatePlace();
	MDOUBLE oldL=VERYSMALL;
	MDOUBLE currL = VERYSMALL;
	tree oldT = _et;
	for (int i=0; i < maxIterations; ++i) {
		if(_unObservableData_p)
			_unObservableData_p->setLforMissingData(_et,&_sp);
		computeUp();
		currL = likelihoodComputation::getTreeLikelihoodFromUp2(_et,_sc,_sp,_cup,_posLike,_weights,_unObservableData_p);
		oldT = _et;
		if (currL < oldL + epsilon) { // need to break
			if (currL<=oldL) {
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
	if(_unObservableData_p)
		_unObservableData_p->setLforMissingData(_et,&_sp);
	currL = likelihoodComputation::getTreeLikelihoodFromUp2(_et,_sc,_sp,_cup,_posLike,_weights, _unObservableData_p);
	//////////////
	//MDOUBLE checkUpLL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et, _sc, _sp, _weights, _plogLforMissingData);
	//LOGnOUT(4, << "checkUpLL = "<<checkUpLL <<" curll = "<<currL<<endl);
	///////////////

	if (currL<=oldL) 
	{
		_et = oldT;
		if(_unObservableData_p)
			_unObservableData_p->setLforMissingData(_et,&_sp);
		return oldL; // keep the old tree, and old likelihood
	} 
	else 
        return currL;
}

void bblEMfixRoot::allocatePlace() {
	_computeCountsV.resize(_et.getNodesNum());//initiateTablesOfCounts
	for (int i=0; i < _computeCountsV.size(); ++i) {
	{
		_computeCountsV[i].resize(_sp.alphabetSize()); 
		for (int letterAtRoot = 0; letterAtRoot < _computeCountsV[0].size(); ++letterAtRoot)
			_computeCountsV[i][letterAtRoot].countTableComponentAllocatePlace(_sp.alphabetSize(),_sp.categories());
			//_computeCountsV[i][letterAtRoot].zero();
		}
	}
	_cup.allocatePlace(_sc.seqLen(),_sp.categories(), _et.getNodesNum(), _sc.alphabetSize());
	_cdown.resize(_sp.categories());
	for (int categor = 0; categor < _sp.categories(); ++categor)
	{
		_cdown[categor].allocatePlace(_sp.alphabetSize(), _et.getNodesNum(), _sc.alphabetSize());
	}	
}

void bblEMfixRoot::bblEM_it(const MDOUBLE tollForPairwiseDist){
	for (int i=0; i < _computeCountsV.size(); ++i) {
		for (int j=0; j < _computeCountsV[0].size(); ++j) {
			_computeCountsV[i][j].zero();
		}
	}
	for (int i=0; i < _sc.seqLen(); ++i) {
		computeDown(i);
		addCounts(i); // computes the counts and adds to the table.
	}
	optimizeBranches(tollForPairwiseDist);
	if(_unObservableData_p)
		_unObservableData_p->setLforMissingData(_et,&_sp);
}

void bblEMfixRoot::optimizeBranches(const MDOUBLE tollForPairwiseDist){
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			fromCountTableComponentToDistancefixRoot from1(_computeCountsV[mynode->id()],_sp,tollForPairwiseDist,mynode->dis2father(),_unObservableData_p);
			from1.computeDistance();
			mynode->setDisToFather(from1.getDistance());
			if(_unObservableData_p)
				_unObservableData_p->setLforMissingData(_et,&_sp);
		}
	}
}

void bblEMfixRoot::computeUp(){
	_pij.fillPij(_et,_sp,0); // 0 is becaues we compute Pij(t) and not its derivations...
	computeUpAlg cupAlg;
	for (int pos=0; pos < _sc.seqLen(); ++pos) {
        for (int categor = 0; categor < _sp.categories(); ++categor) {
			cupAlg.fillComputeUp(_et,_sc,pos,_pij[categor],_cup[pos][categor]);
		}
	}
 }

void bblEMfixRoot::computeDown(const int pos){
	computeDownAlg cdownAlg;
		for (int categor = 0; categor < _sp.categories(); ++categor) {
			cdownAlg.fillComputeDownNonReversible(_et,_sc,pos,_pij[categor],_cdown[categor],_cup[pos][categor]);
		}
 }

void bblEMfixRoot::addCounts(const int pos){
	//MDOUBLE posProb = 
	//	likelihoodComputation::getProbOfPosWhenUpIsFilledGam(pos,_et,_sc,_sp,_cup);
						
	MDOUBLE weig = (_weights ? (*_weights)[pos] : 1.0);
	if (weig == 0) return;
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			addCountsFixedRoot(pos,mynode,_posLike[pos],weig);
		}
	}
}

void bblEMfixRoot::addCountsFixedRoot(const int pos, tree::nodeP mynode, const doubleRep posProb, const MDOUBLE weig){

	computeCounts cc;
	for(int letterAtRoot = 0; letterAtRoot < _sp.alphabetSize(); letterAtRoot++)
	{
		for (int categor =0; categor< _sp.categories(); ++ categor) 
		{
				cc.computeCountsNodeFatherNodeSonHomPos(_sc,
											_pij[categor],
											_sp,
											_cup[pos][categor],
											_cdown[categor][letterAtRoot],
											weig,
											posProb,
											mynode,
											_computeCountsV[mynode->id()][letterAtRoot][categor],
											_sp.ratesProb(categor),
											letterAtRoot);	// letterInFather is used - FixedRoot version
		}
	}
}          
