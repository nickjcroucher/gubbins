// 	$Id: bblEM2USSRV.h 1504 2007-01-15 14:04:44Z osnatz $	
//copy of bblEM of the codon model + changes
#ifndef ___BBL_EM_2_USSRV
#define ___BBL_EM_2_USSRV

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "countTableComponent.h"
#include "computePijComponent.h"
#include "suffStatComponent.h"
#include "ussrvModel.h"
#include "computeUpAlg.h"
#include "computeDownAlg.h"
#include "computeCounts.h"
#include "treeIt.h"
#include "fromCountTableComponentToDistance2USSRV.h"
#include "likelihoodComputation2USSRV.h"
#include "someUtil.h"
#include <vector>
using namespace std;
// @@@@ maybe should inherit from bblEM
class bblEM2USSRV {
public:
	explicit bblEM2USSRV(tree& et,
				const sequenceContainer& sc,
				const sequenceContainer& baseSc,
				const ussrvModel &model,
				const Vdouble * weights = NULL,
				const int maxIterations=50,
				const MDOUBLE epsilon=0.05,
				const MDOUBLE tollForPairwiseDist=0.001);
	MDOUBLE getTreeLikelihood() const {return _treeLikelihood;}

private:
	MDOUBLE compute_bblEM(int maxIterations,
					MDOUBLE epsilon,
					MDOUBLE tollForPairwiseDist);
	void bblEM_it(MDOUBLE tollForPairwiseDist);
	void computeDown(int pos);
	void computeUp();
	void addCounts(int pos);
	void addCounts(int pos, tree::nodeP mynode, doubleRep posProb, MDOUBLE weig);
	void optimizeBranches(MDOUBLE tollForPairwiseDist);
	void allocatePlace();

	MDOUBLE _treeLikelihood;
	tree& _et;
	const sequenceContainer& _sc;
	const sequenceContainer& _baseSc;
	const ussrvModel& _model;
	vector<countTableComponentGam> _computeCountsBaseV; // for each node - a table of rate*alph*alph (see below)
	vector<countTableComponentHom> _computeCountsSsrvV; // for each node - a table of rate*alph*alph (see below)
	computePijGam _pijBase;
	computePijHom _pijSSRV;
	suffStatGlobalGam _cupBase;
	suffStatGlobalHom _cupSSRV; 
	suffStatGlobalGamPos _cdownBase;
	suffStatGlobalHomPos _cdownSSRV;
	const Vdouble * _weights;
	VdoubleRep _posLike;
};

// _computeCountsV is a vector containing for each node a countTableComponentGam.
// countTableComponentGam is a vector containing for each rate category a table of size alphabet*alphabet
// (VVdouble) which should be pre-filled with Pij(x,y,rk) from equation (17) in the EM-BBL theory summary.
// Pij(x,y,rk) represents the probability of observing x and y along a branch ti at position j with rate from 
// category k. 
// For this reason, we need to initialize this class and calculate it again for every position.


#endif // bblEM2USSRV
