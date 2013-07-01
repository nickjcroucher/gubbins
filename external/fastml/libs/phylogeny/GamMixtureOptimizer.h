#ifndef __GAMMIXTURE_OPTIMIZER
#define __GAMMIXTURE_OPTIMIZER
/************************************************************
GamMixtureOptimizer class is used to find the best Gamma mixture parameters.
The parameters to otimized are the alpha and beta of each component and the components probabilities.
The optimizer can choose between several optimization algorithms (EM, ConjugateDerivatives, etc).
The interface to the optimizer is the functions:
1. findBestParam() = given a gammaMixture - finds the best parameters.
2. findBestParamManyStarts() - finds the best parameters but starts from many initial points.
3. SetOptAlg() - choose the optimization algorithm to be used.   
************************************************************/
#include "definitions.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "tree.h"
#include "mixtureDistribution.h"
#include "unObservableData.h"



class GamMixtureOptimizer{
public:
	enum OptimAlg {EM, ONE_DIM, TX_CONJUGATE_DERIVATIVES, NR_CONJUGATE_DERIVATIVES};
public:

	explicit GamMixtureOptimizer(stochasticProcess* cur_sp, const sequenceContainer& sc, const tree& inTree, unObservableData* unObservableData_p = NULL);
	virtual ~GamMixtureOptimizer();

	const stochasticProcess* getSp()  const {return _pSp;}
	const mixtureDistribution* getMixtureDist()  const {return static_cast<mixtureDistribution*>(_pSp->distr());}

	MDOUBLE findBestParamManyStarts(const Vint pointsNum, const Vint iterNum, const vector<OptimAlg> OptAlgs, const Vdouble tols, const Vdouble * pWeights, ofstream* pOutF = NULL);
	//return the logLikelihood. the final distribution is stored in the stochasticProcess
	MDOUBLE findBestParam(const OptimAlg alg, const int maxIterations, const MDOUBLE tol, const Vdouble * pWeights, ofstream* pOutF=NULL);

	void setTolOptSpecific(const MDOUBLE tol) {_tolOptSpecific = tol;}

private:
	MDOUBLE optimizeParam(mixtureDistribution* pInDistribution, const int maxIterations, const OptimAlg alg, const MDOUBLE tol, const Vdouble * pWeights, ofstream* pOutF);


private:	
	stochasticProcess* _pSp;
	const sequenceContainer* _pSc;
	const tree* _pTree;
	unObservableData* _unObservableData_p;

	MDOUBLE _tolOptSpecific; //tolerance specific to the optimization algorithm
};

#endif

