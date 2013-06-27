// $Id: bblEM2codon.cpp 2350 2007-08-20 10:53:51Z adist $
#include "bblEM2codon.h"
#include "likelihoodComputation.h"
#include "likelihoodComputation2Codon.h"
#include "fromCountTableComponentToDistance2Codon.h"
using namespace likelihoodComputation;
using namespace likelihoodComputation2Codon;
#include "computeUpAlg.h"
#include "computeDownAlg.h"
#include "computeCounts.h"
#include "treeIt.h"
#include "errorMsg.h"
#include "logFile.h"
#include <ctime>

bblEM2codon::bblEM2codon(tree& et,
				const sequenceContainer& sc,
				const vector<stochasticProcess>& spVec,
				const distribution *in_distr,
				const Vdouble * weights,
				const int maxIterations,
				const MDOUBLE epsilon,
				const MDOUBLE tollForPairwiseDist) :
  _et(et),_sc(sc),_spVec(spVec),_distr(in_distr->clone()),_weights (weights) {
	
	LOG(5,<<"******BEGIN OF BBL EM*********"<<endl<<endl);
	_treeLikelihood = compute_bblEM(maxIterations,epsilon,tollForPairwiseDist);
	LOG(5,<<"******END OF BBL EM*********"<<endl<<endl);
}

bblEM2codon::~bblEM2codon(){  
 delete _distr;
 }

MDOUBLE bblEM2codon::compute_bblEM(
			const int maxIterations,
			const MDOUBLE epsilon,
			const MDOUBLE tollForPairwiseDist){
	allocatePlace();
	MDOUBLE oldL=VERYSMALL;
	MDOUBLE currL = VERYSMALL;
	tree oldT = _et;
	for (int i=0; i < maxIterations; ++i) {
		
		computeUp();
		//currL = likelihoodComputation::getTreeLikelihoodFromUp2(_et,_sc,_sp,_cup,_posLike,_weights);
		currL = likelihoodComputation2Codon::getTreeLikelihoodFromUp2(_et,_sc,_spVec[0],_cup,_posLike,_distr,_weights);
		//////////////
		if (i!=0)
		  LOG(5,<<"last best L= "<<oldL<<endl);
		LOG(5,<<"current best L= "<<currL<<endl<<endl);
	
		//MDOUBLE checkUpLL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et, _sc, _sp, _weights);
		//cerr << "checkUpLL = "<<checkUpLL <<" curll = "<<currL<<endl;
		///////////////
		
		if (currL < oldL + epsilon) { // need to break
			if (currL<oldL) {
				_et = oldT;
				return oldL; // keep the old tree, and old likelihood
			} else {
                //update the tree and likelihood and return
				return currL;
			}
		}
		oldT = _et;
		bblEM_it(tollForPairwiseDist);
		oldL = currL;
	}
	// in the case were we reached max_iter, we have to recompute the likelihood of the new tree...
	computeUp();
	currL = likelihoodComputation2Codon::getTreeLikelihoodFromUp2(_et,_sc,_spVec[0],_cup,_posLike,_distr,_weights);
	//currL = likelihoodComputation::getTreeLikelihoodFromUp2(_et,_sc,_sp,_cup,_posLike,_weights);
	if (currL<oldL) {
		_et = oldT;
		return oldL; // keep the old tree, and old likelihood
	} 
	else 
        return currL;
}

void bblEM2codon::allocatePlace() {
	_computeCountsV.resize(_et.getNodesNum()); //initiateTablesOfCounts
	for (int i=0; i < _computeCountsV.size(); ++i) {
		_computeCountsV[i].countTableComponentAllocatePlace(_spVec[0].alphabetSize(),_distr->categories());
	}
	_cup.allocatePlace(_sc.seqLen(),_distr->categories(), _et.getNodesNum(), _sc.alphabetSize());
	_cdown.allocatePlace(_distr->categories(), _et.getNodesNum(), _sc.alphabetSize());
}

void bblEM2codon::bblEM_it(const MDOUBLE tollForPairwiseDist){
	int i;
	for (i=0; i < _computeCountsV.size(); ++i) {
		_computeCountsV[i].zero();
	}
	for (i=0; i < _sc.seqLen(); ++i) {
		computeDown(i);
		addCounts(i); // computes the counts and adds to the table.
	}
	optimizeBranches(tollForPairwiseDist);
}

void bblEM2codon::optimizeBranches(const MDOUBLE tollForPairwiseDist){
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			fromCountTableComponentToDistance2Codon from1(_computeCountsV[mynode->id()],_spVec,tollForPairwiseDist,mynode->dis2father());
			from1.computeDistance();
			mynode->setDisToFather(from1.getDistance());
		}
	}
}

void bblEM2codon::computeUp(){
	//_pij.fillPij(_et,_sp,0); // 0 is becaues we compute Pij(t) and not its derivations...
	_pij._V.resize(_spVec.size());
	for (int i=0; i < _spVec.size(); ++i) {
		_pij._V[i].fillPij(_et,_spVec[i]);
	}
	computeUpAlg cupAlg;
	for (int pos=0; pos < _sc.seqLen(); ++pos) {
		for (int categor = 0; categor < _spVec.size(); ++categor) {
			cupAlg.fillComputeUp(_et,_sc,pos,_pij[categor],_cup[pos][categor]);
		}
	}
 }

void bblEM2codon::computeDown(const int pos){
	computeDownAlg cdownAlg;
	for (int categor = 0; categor < _distr->categories(); ++categor) {
		cdownAlg.fillComputeDown(_et,_sc,pos,_pij[categor],_cdown[categor],_cup[pos][categor]);
	}
 }

void bblEM2codon::addCounts(const int pos){
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

void bblEM2codon::addCounts(const int pos, tree::nodeP mynode, const MDOUBLE posProb, const MDOUBLE weig){

	computeCounts cc;
	for (int categor =0; categor< _distr->categories(); ++ categor) {
			cc.computeCountsNodeFatherNodeSonHomPos(_sc,
										_pij[categor],
										_spVec[categor],
										_cup[pos][categor],
										_cdown[categor],
										weig,
										posProb,
										mynode,
										_computeCountsV[mynode->id()][categor],
										_distr->ratesProb(categor));
	}
}          

