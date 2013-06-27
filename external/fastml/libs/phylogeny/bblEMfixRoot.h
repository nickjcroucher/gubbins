// $Id: bblEM.h 4478 2008-07-17 17:09:55Z cohenofi $
#ifndef ___BBL_EM_GL__FIXED_ROOT
#define ___BBL_EM_GL__FIXED_ROOT

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "countTableComponent.h"
#include "computePijComponent.h"
#include "suffStatComponent.h"
#include "gainLossAlphabet.h"
#include "unObservableData.h"
#include <vector>

using namespace std;

class bblEMfixRoot {
public:
	explicit bblEMfixRoot(tree& et,
				const sequenceContainer& sc,
				const stochasticProcess& sp,
				const Vdouble * weights = NULL,
				const int maxIterations=50,
				const MDOUBLE epsilon=0.05,
				const MDOUBLE tollForPairwiseDist=0.001,
				unObservableData*  _unObservableData_p=NULL);
	MDOUBLE getTreeLikelihood() const {return _treeLikelihood;}

private:
	MDOUBLE compute_bblEM(const int maxIterations,
					const MDOUBLE epsilon,
					const MDOUBLE tollForPairwiseDist);
	void bblEM_it(const MDOUBLE tollForPairwiseDist);
	void computeDown(const int pos);
	void computeUp();
	void addCounts(const int pos);
	void addCountsFixedRoot(const int pos, tree::nodeP mynode, const doubleRep posProb, const MDOUBLE weig);

	void optimizeBranches(const MDOUBLE tollForPairwiseDist);
	void allocatePlace();



	MDOUBLE _treeLikelihood;
	tree& _et;
	const sequenceContainer& _sc;
	const stochasticProcess& _sp;
	//vector<countTableComponentGam> _computeCountsV; // for each node - a table of rate*alph*alph
	vector< vector< countTableComponentGam > > _computeCountsV; // _computeCountsV[node][letterAtRoot][rate][alph][alph]
	computePijGam _pij;
	suffStatGlobalGam _cup;						//_cup[pos][categ][nodeid][letter][prob]			
	//suffStatGlobalGamPos _cdown;				// foreach pos: computeDown(pos);	addCounts(pos);  
	vector<suffStatGlobalGamPos> _cdown;		//_cdown[categ][letter@root][nodeid][letter][prob] - since fillComputeDownNonReversible uses this assumption
	const Vdouble * _weights;
	VdoubleRep _posLike;
	unObservableData*  _unObservableData_p;
};

#endif
