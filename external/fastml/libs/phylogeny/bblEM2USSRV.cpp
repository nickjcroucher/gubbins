// 	$Id: bblEM2USSRV.cpp 1944 2007-04-18 12:41:14Z osnatz $	
#include "bblEM2USSRV.h"

bblEM2USSRV::bblEM2USSRV(tree& et,
				const sequenceContainer& sc,
				const sequenceContainer& baseSc,
				const ussrvModel& model,
				const Vdouble * weights,
				int maxIterations,
				MDOUBLE epsilon,
				MDOUBLE tollForPairwiseDist) :
_et(et),_sc(sc),_baseSc(baseSc),_model(model),_weights (weights)
{	
	LOG(5,<<"******BBL EM USSRV*********"<<endl<<endl);
	_treeLikelihood = compute_bblEM(maxIterations,epsilon,tollForPairwiseDist);
}

// @@@@ Need to check if we can make it more efficient
MDOUBLE bblEM2USSRV::compute_bblEM(
			int maxIterations,
			MDOUBLE epsilon,
			MDOUBLE tollForPairwiseDist){
	
	allocatePlace();
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE currL = VERYSMALL;
	tree oldT = _et;
	for (int i=0; i < maxIterations; ++i) {
		computeUp();
		// Calculate the likelihood and fill the _posLike
		currL = likelihoodComputation2USSRV::getTreeLikelihoodFromUp2(_et,
				_sc,_baseSc,_model,_cupBase,_cupSSRV,_posLike,_weights);
		//////////////
		LOGDO(5,printTime(myLog::LogFile()));
		LOG(5,<<"iteration no "<<i << " in BBL "<<endl);
		LOG(5,<<"old best  L= "<<oldL<<endl);
		LOG(5,<<"current best  L= "<<currL<<endl);
	

		if (currL < oldL + epsilon) { // need to break
			if (currL<oldL) {
				cout<<"******** PROBLEMS IN BBL USSRV*********"<<endl;
				LOG(5,<<"old best  L= "<<oldL<<endl);
				LOG(5,<<"current best  L= "<<currL<<endl);
				_et = oldT;
				return oldL; // keep the old tree, and old likelihood
			} else {
                //update the tree and likelihood and return
				LOG(5,<<"old best  L= "<<oldL<<endl);
				LOG(5,<<"current best  L= "<<currL<<endl);
				return currL;
			}
		}
		oldT = _et;
		bblEM_it(tollForPairwiseDist);
		oldL = currL;
	}
	// in the case were we reached max_iter, we have to recompute the likelihood of the new tree...
	computeUp();
	currL = likelihoodComputation2USSRV::getTreeLikelihoodFromUp2(_et,
			_sc,_baseSc,_model,_cupBase,_cupSSRV,_posLike,_weights);
	if (currL<oldL) {
		_et = oldT;
		return oldL; // keep the old tree, and old likelihood
	} 
	else 
        return currL;
}


void bblEM2USSRV::allocatePlace() {
	_computeCountsBaseV.resize(_et.getNodesNum()); //initiateTablesOfCounts
	_computeCountsSsrvV.resize(_et.getNodesNum()); //initiateTablesOfCounts
	
	for (int i=0; i < _computeCountsBaseV.size(); ++i) {
		_computeCountsBaseV[i].countTableComponentAllocatePlace(_model.getBaseModel().alphabetSize(),_model.noOfCategor());
		_computeCountsSsrvV[i].countTableComponentAllocatePlace(_model.getSSRVmodel().alphabetSize());
	}
	_cupBase.allocatePlace(_baseSc.seqLen(),_model.noOfCategor(), _et.getNodesNum(), _baseSc.alphabetSize());
	_cupSSRV.allocatePlace(_sc.seqLen(), _et.getNodesNum(), _sc.alphabetSize());

	_cdownBase.allocatePlace(_model.noOfCategor(), _et.getNodesNum(), _baseSc.alphabetSize());
	_cdownSSRV.allocatePlace( _et.getNodesNum(), _sc.alphabetSize());

}

void bblEM2USSRV::bblEM_it(MDOUBLE tollForPairwiseDist){
	for (int i=0; i < _computeCountsBaseV.size(); ++i) {
		_computeCountsBaseV[i].zero();
		_computeCountsSsrvV[i].zero();
	}
	for (int i=0; i < _sc.seqLen(); ++i) {
		computeDown(i);
		addCounts(i); // computes the counts and adds to the table.
	}
	optimizeBranches(tollForPairwiseDist);
}

// @@@@ need to print the tree
void bblEM2USSRV::optimizeBranches(MDOUBLE tollForPairwiseDist){
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			fromCountTableComponentToDistance2USSRV 
				from1(_computeCountsBaseV[mynode->id()],_computeCountsSsrvV[mynode->id()],_model,tollForPairwiseDist,mynode->dis2father());
			from1.computeDistance();
			mynode->setDisToFather(from1.getDistance());
		}
	}
}

void bblEM2USSRV::computeUp(){
	_pijBase.fillPij(_et,_model.getBaseModel(),0); // 0 is becaues we compute Pij(t) and not its derivations...
	_pijSSRV.fillPij(_et,_model.getSSRVmodel(),0);
	
	computeUpAlg cupAlg;
	for (int pos=0; pos < _sc.seqLen(); ++pos) {
		// compute up for the base model
		for (int categor = 0; categor < _model.noOfCategor(); ++categor) {
			cupAlg.fillComputeUp(_et,_baseSc,pos,_pijBase[categor],_cupBase[pos][categor]);
		}
		// compute up for the ssrv model
		cupAlg.fillComputeUp(_et,_sc,pos,_pijSSRV,_cupSSRV[pos]);
	}
}

void bblEM2USSRV::computeDown(int pos){
	computeDownAlg cdownAlg;
	// compute down for the base model
	for (int categor = 0; categor < _model.noOfCategor(); ++categor) {
		cdownAlg.fillComputeDown(_et,_baseSc,pos,_pijBase[categor],_cdownBase[categor],_cupBase[pos][categor]);		
	}
	// compute down for the ssrv model
	cdownAlg.fillComputeDown(_et,_sc,pos,_pijSSRV,_cdownSSRV,_cupSSRV[pos]);		
}

void bblEM2USSRV::addCounts(int pos){
						
	MDOUBLE weig = (_weights ? (*_weights)[pos] : 1.0);
	if (weig == 0) return;
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			addCounts(pos,mynode,_posLike[pos],weig);
		}
	}
}

void bblEM2USSRV::addCounts(int pos, tree::nodeP mynode, doubleRep posProb, MDOUBLE weig){

	computeCounts cc;
	int categor;
	// base Model
	for (categor =0; categor< _model.noOfCategor(); ++categor) {
			cc.computeCountsNodeFatherNodeSonHomPos(_baseSc, 
										_pijBase[categor],
										_model.getBaseModel(),
										_cupBase[pos][categor],
										_cdownBase[categor],
										weig,
										posProb,
										mynode,
										_computeCountsBaseV[mynode->id()][categor],
										_model.getCategorProb(categor)*(1-_model.getF()));
	
	}
	// SSRV model
	cc.computeCountsNodeFatherNodeSonHomPos(_sc, 
										_pijSSRV,
										_model.getSSRVmodel(),
										_cupSSRV[pos],
										_cdownSSRV,
										weig,
										posProb,
										mynode,
										_computeCountsSsrvV[mynode->id()],
										_model.getF());
}          



