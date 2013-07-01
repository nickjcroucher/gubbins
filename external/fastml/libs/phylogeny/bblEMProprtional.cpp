// $Id: bblEMProprtional.cpp 962 2006-11-07 15:13:34Z privmane $
#include "bblEM.h"
#include "bblEMProportional.h"
#include "likelihoodComputation.h"
using namespace likelihoodComputation;
#include "computeUpAlg.h"
#include "computeDownAlg.h"
#include "computeCounts.h"
#include "treeIt.h"
#include "fromCountTableComponentToDistance.h"
#include <ctime>//#define VERBOS
#include "fromCountTableComponentToDistanceProp.h"

bblEMProportional::bblEMProportional(tree& et,
									const vector<sequenceContainer>& sc,
									const vector<stochasticProcess>& sp,
									const vector<Vdouble *> * weights,
									const int maxIterations,
									const MDOUBLE epsilon,
									const MDOUBLE tollForPairwiseDist):

_et(et),_sc(sc),_sp(sp),_weights (weights) {
	_numberOfGenes = _sc.size();
	assert(_sp.size() == _sc.size());
	_treeLikelihood = compute_bblEMProp(maxIterations,epsilon,tollForPairwiseDist);
}

MDOUBLE bblEMProportional::compute_bblEMProp(
			const int maxIterations,
			const MDOUBLE epsilon,
			const MDOUBLE tollForPairwiseDist){
	allocatePlaceProp();
	MDOUBLE oldL=VERYSMALL;
	MDOUBLE currL = VERYSMALL;
	for (int i=0; i < maxIterations; ++i) {
		computeUpProp();
		currL = 0;
		for (int geneN=0; geneN < _numberOfGenes; ++geneN) {
			currL += likelihoodComputation::getTreeLikelihoodFromUp2(_et,_sc[geneN],_sp[geneN],_cup[geneN],_posLike[geneN],(_weights?(*_weights)[geneN]:NULL));
		}
		tree oldT = _et;
		if (currL < oldL + epsilon) { // need to break
			if (currL<oldL) {
				_et = oldT;
				return oldL; // keep the old tree, and old likelihood
			} else {
                //update the tree and likelihood and return
				return currL;
			}
		}
		bblEM_itProp(tollForPairwiseDist);
		oldL = currL;
	}
	return currL;
}

void bblEMProportional::allocatePlaceProp() {
	_computeCountsV.resize(_numberOfGenes);
	_cup.resize(_numberOfGenes);
	_cdown.resize(_numberOfGenes);
	_pij.resize(_numberOfGenes);
	_posLike.resize(_numberOfGenes);
	for (int geneN=0; geneN < _numberOfGenes; ++geneN) {
		_computeCountsV[geneN].resize(_et.getNodesNum()); //initiateTablesOfCounts
		for (int i=0; i < _computeCountsV[geneN].size(); ++i) {
			_computeCountsV[geneN][i].countTableComponentAllocatePlace(_sp[geneN].alphabetSize(),_sp[geneN].categories());
		}
		_cup[geneN].allocatePlace(_sc[geneN].seqLen(),_sp[geneN].categories(), _et.getNodesNum(), _sc[geneN].alphabetSize());
		_cdown[geneN].allocatePlace(_sp[geneN].categories(), _et.getNodesNum(), _sc[geneN].alphabetSize());
	}
}

void bblEMProportional::computeUpProp(){
	for (int geneN=0; geneN < _numberOfGenes; ++geneN) {
		_pij[geneN].fillPij(_et,_sp[geneN],0); // 0 is becaues we compute Pij(t) and not its derivations...
		computeUpAlg cupAlg;
		for (int pos=0; pos < _sc[geneN].seqLen(); ++pos) {
			for (int categor = 0; categor < _sp[geneN].categories(); ++categor) {
				cupAlg.fillComputeUp(_et,_sc[geneN],pos,_pij[geneN][categor],_cup[geneN][pos][categor]);
			}
		}
	}
 }

void bblEMProportional::bblEM_itProp(const MDOUBLE tollForPairwiseDist){
	for (int geneN=0; geneN < _numberOfGenes; ++geneN) {
		for (int i=0; i < _computeCountsV.size(); ++i) {
		 	_computeCountsV[geneN][i].zero();
		}
		for (int i=0; i < _sc[geneN].seqLen(); ++i) {
			computeDownProp(geneN,i);
			addCountsProp(geneN,i); // computes the counts and adds to the table.
		}
	}
	optimizeBranchesProp(tollForPairwiseDist);
}

void bblEMProportional::computeDownProp(const int gene, const int pos){
	computeDownAlg cdownAlg;
	for (int categor = 0; categor < _sp[gene].categories(); ++categor) {
		cdownAlg.fillComputeDown(_et,_sc[gene],pos,_pij[gene][categor],_cdown[gene][categor],_cup[gene][pos][categor]);
	}
}

void bblEMProportional::addCountsProp(const int gene, const int pos){
	vector<MDOUBLE> * weightsOfGene = (_weights?(*_weights)[gene]:NULL);					
	MDOUBLE weig = (weightsOfGene ? (*weightsOfGene)[pos] : 1.0);
	if (weig == 0) return;
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			addCountsProp(gene,pos,mynode,_posLike[gene][pos],weig);
		}
	}
}

void bblEMProportional::addCountsProp(const int gene,const int pos, tree::nodeP mynode, const doubleRep posProb, const MDOUBLE weig){
	computeCounts cc;
	for (int categor =0; categor< _sp[gene].categories(); ++ categor) {
			cc.computeCountsNodeFatherNodeSonHomPos(_sc[gene],
										_pij[gene][categor],
										_sp[gene],
										_cup[gene][pos][categor],
										_cdown[gene][categor],
										weig,
										posProb,
										mynode,
										_computeCountsV[gene][mynode->id()][categor],
										_sp[gene].ratesProb(categor));
	}
}

void bblEMProportional::optimizeBranchesProp(const MDOUBLE tollForPairwiseDist){
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			fromCountTableComponentToDistanceProp from1(_computeCountsV[mynode->id()],_sp,tollForPairwiseDist,mynode->dis2father());
			from1.computeDistance();
			mynode->setDisToFather(from1.getDistance());
		}
	}
}
