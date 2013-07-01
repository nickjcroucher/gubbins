#ifndef ___OPT_GAMMA_MIXTURE_EM
#define ___OPT_GAMMA_MIXTURE_EM
/************************************************************
optGammaMixtureEM class is used to maximize the gammaMixture parameters.
The parameters to otimized are the alpha and beta of each component and the components probabilities.
In each iteration: 
(1) The sufficient statistics are calculated.
(2) Based on these statistics the parameters are optimized.
the procedure stops when no improvment in the tree likelihood is achieved
************************************************************/
#include "definitions.h"
#include "suffStatGammaMixture.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "tree.h"
#include "gammaUtilities.h"

#include <cmath>

class optGammaMixtureEM{

public:
	explicit optGammaMixtureEM(const stochasticProcess& cur_sp, const sequenceContainer& sc, const tree& inTree);
	virtual ~optGammaMixtureEM();

	//return the logLikelihood. the final distribution is stored in the stochasticProcess
	MDOUBLE optimizeParam(mixtureDistribution* pInDistribution, const int maxIterations, const MDOUBLE epsilon, const MDOUBLE epsilomQopt, ofstream* pOutF);

	const stochasticProcess* getSp()  const {return _pSp;}

	MDOUBLE findBestParamManyStarts(const int startPointsNum, const int bestStartsNum, const int startIter, const int maxIterations, const MDOUBLE epsilon, const MDOUBLE epsilomQopt, ofstream* pOutF = NULL);


	void maximizeGammaParam(stochasticProcess* pNewSp, MDOUBLE accuracy);
	void maximizeGammaParam(const suffStatGammaMixture & stats, stochasticProcess* pNewSp, MDOUBLE accuracy);
private:
	void printIter(const stochasticProcess& pInSp, const int it, const MDOUBLE curL); 


	MDOUBLE findBestAlpha(const suffStatGammaMixture& stats, const int compNum, const MDOUBLE accuracy, const MDOUBLE upperBoundAlpha) const;

	void checkEntropy(stochasticProcess & oldSp, stochasticProcess & inSp);


private:	
	stochasticProcess* _pSp;
	const sequenceContainer* _pSc;
	const tree* _pTree;
};




class C_evalAlphaEM{
public:
  explicit C_evalAlphaEM(const suffStatGammaMixture& stats, const int compNum)
			:_compNum(compNum) {_pStats = &stats;}

public:
	MDOUBLE operator() (const MDOUBLE x) 
	{
		MDOUBLE Ak = _pStats->getAk(_compNum);
		MDOUBLE Bk = _pStats->getBk(_compNum);
		MDOUBLE Mk = _pStats->getMk(_compNum);

		MDOUBLE res = log(x) - gammaDerivative(x) + log(Mk) - log(Ak) + (Bk / Mk);
		//cerr<<"+++++++ x = "<<x<<" Ak = "<<Ak<<" Bk = "<<Bk<<" Mk = "<<Mk<<" RES = "<<res<<endl;
// when x is beta (checking)
//		MDOUBLE res = Mk * log(x) - Mk * diGamma(Ak * x / Mk) + Bk;
		return res;
	}

private:
	MDOUBLE diGammaPlus(MDOUBLE x) const 
	{
		MDOUBLE res1 = log(x) + (1/(2*x)) - (1/(12*x*x)) + (1/(120*pow(x, 4))) - (1/(252*pow(x, 6))); 
		MDOUBLE res = log(x) + (0.5/x) - (0.083333333333333333333333333333333/(x*x)) + (0.0083333333333333333333333333333333/(x*x*x*x)) - (0.003968253968253968253968253968254/(pow(x, 6))); 
		return res;
	}
	MDOUBLE diGamma(MDOUBLE x) const 
	{
		//if x<1: use the identity digamma(Z) = digamma(z+1)- (1/z) see http://mathworld.wolfram.com/DigammaFunction.html 
		if (x < 1) 
			return (diGamma(x+1) - (1.0 / x)); 
		MDOUBLE res = log(x) - (1/(2*x)) - (1/(12*x*x)) + (1/(120*pow(x, 4))) - (1/(252*pow(x, 6)));
		//using more terms in the series expansion:
		MDOUBLE debugRes = log(x) - (1/(2*x)) - (1/(12*x*x)) + (1/(120*pow(x, 4))) - (1/(252*pow(x, 6))) + (1/(240*pow(x, 8))) - (1/(132*pow(x, 10))); 
		return res;
	}

	MDOUBLE gammaDerivative(MDOUBLE x) const 
	{
		//MDOUBLE resCheck = (gammln(x+0.001) - gammln(x)) /0.001;
		MDOUBLE res = diGamma(x);
		return res;
	}
private:
	const suffStatGammaMixture* _pStats;
	const int _compNum;
};
#endif

