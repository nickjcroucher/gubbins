#ifndef ___BEST_ALPHA_AND_K
#define ___BEST_ALPHA_AND_K

#include "definitions.h"
#include "tree.h"
#include "likelihoodComputation.h"
#include "likelihoodComputation2Codon.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "generalGammaDistribution.h"
#include "logFile.h"
#include "wYangModel.h"
#include "bblEM2codon.h"
#include "computeUpAlg.h"
#include "numRec.h"



//evaluate best parameters
class optimizeSelectonParameters {
public:
	explicit optimizeSelectonParameters(tree& et, 
					   const sequenceContainer& sc,
					   vector<stochasticProcess>& spVec,
					   distribution * distr,
					   bool bblFlag = true,
						bool isGamma = true, bool isBetaProbSet=false,bool isOmegaSet = false,
						bool isKappaSet=false, bool isAlphaSet=false, bool isBetaSet=false,
					   const MDOUBLE upperBoundOnAlpha = 3.0, // changed from 20, Adi S. 2/7/07
				       const MDOUBLE upperBoundOnBeta = 3.0, // changed from 20, Adi S. 2/7/07
					   const MDOUBLE epsilonAlphaOptimization= 0.01,
					   const MDOUBLE epsilonKOptimization=0.01,
					   const MDOUBLE epsilonLikelihoodImprovment= 0.1,
					   const int maxBBLIterations=20,
					   const int maxTotalIterations=20);
	const MDOUBLE getBestAlpha() const{return _bestAlpha;}
	const MDOUBLE getBestBeta() const{return _bestBeta;}
	const MDOUBLE getBestL() const {return _bestL;}
	const MDOUBLE getBestK() const {return _bestK;}
	const MDOUBLE getBestOmega() const {return _bestOmega;}
	const MDOUBLE getBestBetaProb() const {return _bestBetaProb;}
private:
	MDOUBLE _bestAlpha;
	MDOUBLE _bestL;
	MDOUBLE _bestK;
	MDOUBLE _bestBeta;
	MDOUBLE _bestOmega;
	MDOUBLE _bestBetaProb;
};


//The functor to eval likelihood given a change in a parameters
class evalParam{
public:
	explicit evalParam(const tree& et,
				const sequenceContainer& sc,
				vector<stochasticProcess> spVec,
				int alphaOrKs,
				const distribution * in_distr,
				bool isGamma)
    : _et(et),_sc(sc),_spVec(spVec),_alphaOrKs(alphaOrKs),_isGamma(isGamma){_distr=in_distr->clone();};
	MDOUBLE operator()(MDOUBLE param);
	
    virtual ~evalParam();
	evalParam(const evalParam &other);
	void updateAlpha(MDOUBLE param);
	void updateK(MDOUBLE param);
	void updateBeta(MDOUBLE param);
	void updateOmega(MDOUBLE param);
	void updateBetaProb(MDOUBLE param);
private:
	const tree& _et;
	const sequenceContainer& _sc;

	vector<stochasticProcess> _spVec;
	int _alphaOrKs; //flag to eval different parameters (alpha,beta or ks)
	distribution *_distr;
	bool _isGamma; //gamma = true/ beta=false

};

#endif


