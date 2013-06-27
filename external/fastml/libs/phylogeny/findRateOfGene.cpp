// $Id: findRateOfGene.cpp 962 2006-11-07 15:13:34Z privmane $

#include "definitions.h"
#include "findRateOfGene.h"
#include "computeUpAlg.h"

//#define VERBOS

class findRateOfGene{
public:
  explicit findRateOfGene(const tree &t,
						 const sequenceContainer& sc,
						 stochasticProcess& sp,
						 const Vdouble * weights):  _t(t), _sc(sc),
								_sp(sp),_weights(weights){};
private:
	const tree& _t;
	const sequenceContainer& _sc;
	stochasticProcess& _sp;
	const Vdouble * _weights;
public:
	MDOUBLE operator() (const MDOUBLE fac) {
#ifdef VERBOS
	LOG(5,<<"factor = "<<fac<<endl);
#endif
		_sp.setGlobalRate(fac);
		MDOUBLE tmp =   likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_t,_sc,_sp,_weights);
#ifdef VERBOS
	LOG(5,<<"likelihood = "<<tmp<<endl);
#endif
		return -tmp;
	}
};

MDOUBLE findTheBestFactorFor(const tree &t,
						 const sequenceContainer& sc,
						 stochasticProcess& sp,
						 const Vdouble * weights,
						  MDOUBLE & logLresults) {
#ifdef VERBOS
  LOG(5,<<"xxx   in funtion findTheNestFactorFor xxxxxxxxx"<<endl);
  LOG(5,<<"xxx   b4 optimization xxxxxxxxx"<<endl);
  MDOUBLE myL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(t,sc,sp);
  LOG(5,<<" likelihod is: "<<myL<<endl);
  LOG(5,<<" global rate is: "<<sp.getGlobalRate()<<endl);
  LOG(5,<<"\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
#endif

  const MDOUBLE ax=0,bx=1.0,cx=4.0,tol=0.01f;
  MDOUBLE res=-1.0;
  logLresults =-brent(ax,bx,cx,
		  findRateOfGene(t,sc,sp,weights),
		  tol,
		  &res);
#ifdef VERBOS
  LOG(5,<<"rate of gene = "<<res<<endl);
  LOG(5,<<"xxx   in funtion findTheNestFactorFor xxxxxxxxx"<<endl);
  LOG(5,<<"xxx   after optimization xxxxxxxxx"<<endl);
  myL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(t,sc,sp);
  LOG(5,<<" likelihod is: "<<myL<<"\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
#endif
  sp.setGlobalRate(res);
  return res;}

void makeAverageRateEqOne(tree& et,vector<stochasticProcess> & spVec){
	MDOUBLE sumGlobalRates=0.0;
	for (int k=0; k < spVec.size(); ++k) {
		sumGlobalRates+=spVec[k].getGlobalRate();
	}
	for (int j=0; j < spVec.size(); ++j) {
		MDOUBLE newGlobalRate = spVec[j].getGlobalRate();
		newGlobalRate*=(spVec.size()/sumGlobalRates);
		spVec[j].setGlobalRate(newGlobalRate);

	}
	et.multipleAllBranchesByFactor(sumGlobalRates/spVec.size());
}




