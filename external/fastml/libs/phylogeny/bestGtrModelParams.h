// $Id: bestGtrModelparams.h 2008-28-04 15:13:34Z nimrod $

#ifndef ___BEST_GTRMODEL_PARAMS
#define ___BEST_GTRMODEL_PARAMS

#include "definitions.h"

#include "likelihoodComputation.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "gammaDistribution.h"
#include "generalGammaDistribution.h"
#include "tree.h"
#include "gtrModel.h"

typedef enum
  {
    Invalid = 0,
	a2c,
	a2g,
	a2t,
	c2g,
	c2t,
	g2t,
  }GTRParam;

#define maxBBLIt 10
#define epsilonLoglikeForBBL 0.01
#define inAlpha 1.5
#define epsilonLoglikeForAlphaOptimization 0.01
#define upperBoundForAlpha 5.0

class bestGtrModel {
public:
	explicit bestGtrModel(tree& et, // find best Gtr Model Params
										const sequenceContainer& sc,
										stochasticProcess& sp,
										const Vdouble * weights,
										const int maxTotalIterations = 5,
										const MDOUBLE epsilonLikelihoodImprovment = 0.05,
										const MDOUBLE epsilonLoglikelihoodForGTRParam = 0.01,
										const MDOUBLE upperBoundGTRParam = 5.0,
										const bool optimizeTree = true,
                                        const bool optimizeAlpha = true);
	MDOUBLE getBesta2c() {return _best_a2c;}
	MDOUBLE getBesta2g() {return _best_a2g;}
	MDOUBLE getBesta2t() {return _best_a2t;}
	MDOUBLE getBestc2g() {return _best_c2g;}
	MDOUBLE getBestc2t() {return _best_c2t;}
	MDOUBLE getBestg2t() {return _best_g2t;}
	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _best_a2c;
	MDOUBLE _best_a2g;
	MDOUBLE _best_a2t;
	MDOUBLE _best_c2g;
	MDOUBLE _best_c2t;
	MDOUBLE _best_g2t;
	MDOUBLE _bestAlpha;
	MDOUBLE _bestL;
};

class C_evalGTRParam{
public:
  C_evalGTRParam(	const GTRParam param,
					const tree& et,
					const sequenceContainer& sc,
					stochasticProcess& sp,
					const Vdouble * weights = NULL)
		:_param(param), _et(et),_sc(sc),_weights(weights),_sp(sp){};
private:
	const GTRParam _param;
	const tree& _et;
	const sequenceContainer& _sc;
	const Vdouble * _weights;
	stochasticProcess& _sp;
public:
	MDOUBLE operator() (MDOUBLE paramVal) {
		switch (_param){
			case a2c:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_a2c(paramVal);
				break;
			case a2g:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_a2g(paramVal);
				break;	
			case a2t:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_a2t(paramVal);
				break;
			case c2g:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_c2g(paramVal);
				break;
			case c2t:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_c2t(paramVal);
				break;
			case g2t:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_g2t(paramVal);
				break;
			default:
				errorMsg::reportError("Missing GTR parameter in C_evalGTRParam::operator ()");
				break;
		}
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_sp,_weights);
		LOG(5,<<" with " + int2string(_param) + " = "<<paramVal<<" logL = "<<res<<endl);
		return -res;
	}
};

#endif


