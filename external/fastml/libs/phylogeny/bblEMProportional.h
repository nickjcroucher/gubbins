// $Id: bblEMProportional.h 962 2006-11-07 15:13:34Z privmane $
#ifndef ___BBL_EM_PROPORTIONAL_H
#define ___BBL_EM_PROPORTIONAL_H

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"

#include <vector>
using namespace std;


class bblEMProportional {
public:
	explicit bblEMProportional(tree& et,
									const vector<sequenceContainer>& sc,
									const vector<stochasticProcess>& sp,
									const vector<Vdouble *> * weights = NULL,
									const int maxIterations=50,
									const MDOUBLE epsilon=0.05,
									const MDOUBLE tollForPairwiseDist=0.0001);
	MDOUBLE getTreeLikelihood() const {return _treeLikelihood;}

private:
	MDOUBLE compute_bblEMProp(const int maxIterations,const MDOUBLE epsilon,const MDOUBLE tollForPairwiseDist);
	void allocatePlaceProp();
	void computeUpProp();
	void bblEM_itProp(const MDOUBLE tollForPairwiseDist);
	void computeDownProp(const int gene, const int pos);
	void addCountsProp(const int gene, const int pos);
	void addCountsProp(const int gene,const int pos, tree::nodeP mynode, const doubleRep posProb, const MDOUBLE weig);
	void optimizeBranchesProp(const MDOUBLE tollForPairwiseDist);

	MDOUBLE _treeLikelihood;
	tree& _et;
	const vector<sequenceContainer>& _sc;
	const vector<stochasticProcess>& _sp;
	const vector<Vdouble *> * _weights;
	int _numberOfGenes;
	vector<	vector<countTableComponentGam> > _computeCountsV; // for each gene, for each node - a table of rate*alph*alph
	vector<suffStatGlobalGam> _cup;
	vector<suffStatGlobalGamPos> _cdown;
	vector<computePijGam> _pij;
	VVdoubleRep _posLike;


};

#endif
