#ifndef ___OPT_GAMMA_MIXTURE_LS
#define ___OPT_GAMMA_MIXTURE_LS
/************************************************************
optGammaMixtureLS class is used to maximize the gammaMixture parameters via a line search maximization.
The parameters to otimized are the alpha and beta of each component and the components probabilities.
In each iteration: 
optimized all parameters iteratively
The procedure stops when no improvment in the tree likelihood is achieved
************************************************************/
#include "definitions.h"
#include "suffStatGammaMixture.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "tree.h"
#include "gammaUtilities.h"
#include "likelihoodComputation.h"
#include "unObservableData.h"



#include <cmath>

class optGammaMixtureLS{
public:
	enum optAlg {ONE_DIM/*, POWELL, CONJUGATE_DERIVATIVES*/};

public:
	explicit optGammaMixtureLS(stochasticProcess* pSp, const sequenceContainer& sc, const tree& inTree, MDOUBLE upperBoundAlpha =15.0, MDOUBLE upperBoundBeta =15.0, unObservableData* unObservableData_p=NULL);
	virtual ~optGammaMixtureLS();

	//return the logLikelihood. the final distribution is stored in the stochasticProcess
	MDOUBLE optimizeParam(const int maxIterations, const MDOUBLE tol, const Vdouble * pWeights, optAlg optType);
	MDOUBLE optimizeParam(mixtureDistribution * pMixture, const int maxIterations, const MDOUBLE tol, const Vdouble * pWeights, optAlg optType);
	
	
private:
	void printIter(const mixtureDistribution * pMixture, const int it, const MDOUBLE curL); 
	
	MDOUBLE optimizeParamOneDim(mixtureDistribution * pMixture, const int maxIterations, const MDOUBLE tol, const Vdouble * pWeights);
	//MDOUBLE optimizeParamPowell(mixtureDistribution * pMixture, const int maxIterations, const MDOUBLE tol, const Vdouble * pWeights, ofstream* pOutF=NULL);
	//MDOUBLE optimizeParamConjugateDeriv(mixtureDistribution *pMixture, 
	//			const int maxIterations, const MDOUBLE tol, const Vdouble *pWeights, ofstream* pOutF);

	//MDOUBLE optimizeParam1CompPowel(mixtureDistribution * pMixture, const int maxIterations, const MDOUBLE tol, const Vdouble * pWeights, ofstream* pOutF=NULL);
	//MDOUBLE optimizeParamManyCompPowel(mixtureDistribution * pMixture, const int maxIterations, const MDOUBLE tol, const Vdouble * pWeights, ofstream* pOutF=NULL);

private:	
	stochasticProcess* _pSp;
	const sequenceContainer* _pSc;
	const tree* _pTree;
	unObservableData* _unObservableData_p;

	MDOUBLE _upperBoundAlpha;
	MDOUBLE _upperBoundBeta; 
};

//line search classes for brent
class C_evalAlphaMixture{
public:
  C_evalAlphaMixture(const tree& et,
					const sequenceContainer& sc,
					stochasticProcess* pSp,
					const int componetNumber,
					const Vdouble * weights = NULL,
					unObservableData* unObservableData_p=NULL)
		: _et(et),_sc(sc),_weights(weights),_pSp(pSp), _compNum(componetNumber)
  {
	  if(unObservableData_p)
		  _unObservableData_p = unObservableData_p->clone();
	  else
		  _unObservableData_p = NULL;
  };
  virtual ~C_evalAlphaMixture(){
	  if(_unObservableData_p)	delete _unObservableData_p;
  }

private:
	const tree& _et;
	const sequenceContainer& _sc;
	const Vdouble * _weights;
	unObservableData* _unObservableData_p;
	stochasticProcess* _pSp;
	const int _compNum;
public:
	MDOUBLE operator() (MDOUBLE alpha) {
		if (_pSp->categories() == 1) {
			errorMsg::reportError(" one category when trying to optimize alpha");
		}
		mixtureDistribution * pMixture = static_cast<mixtureDistribution*>(_pSp->distr());
		pMixture->setAlpha(alpha, _compNum);
		if(_unObservableData_p){
			_unObservableData_p->setLforMissingData(_et,_pSp);
		}
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,*_pSp,_weights,_unObservableData_p);
#ifdef VERBOS
		cerr<<"Component = "<<_compNum<<" with alpha = "<<alpha<<" logL = "<<res<<endl;
#endif
		return -res;
	}
};

class C_evalBetaMixture{
public:
  C_evalBetaMixture(const tree& et,
					const sequenceContainer& sc,
					stochasticProcess* pSp,
					const int componetNumber,
					const Vdouble * weights = NULL,
					unObservableData* unObservableData_p=NULL)
		: _et(et),_sc(sc),_weights(weights),_pSp(pSp), _compNum(componetNumber)
  {
	  if(unObservableData_p)
		  _unObservableData_p = unObservableData_p->clone();
	  else
		  _unObservableData_p = NULL;
  };
  virtual ~C_evalBetaMixture(){
	  if(_unObservableData_p)	delete _unObservableData_p;
  }

private:
	const tree& _et;
	const sequenceContainer& _sc;
	const Vdouble * _weights;
	unObservableData* _unObservableData_p;
	stochasticProcess* _pSp;
	const int _compNum;
public:
	MDOUBLE operator() (MDOUBLE beta) {
		if (_pSp->categories() == 1) {
			errorMsg::reportError(" one category when trying to optimize beta");
		}
		mixtureDistribution * pMixture = static_cast<mixtureDistribution*>(_pSp->distr());
		pMixture->setBeta(beta, _compNum);
		if(_unObservableData_p){
			_unObservableData_p->setLforMissingData(_et,_pSp);
		}
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,*_pSp,_weights,_unObservableData_p);
#ifdef VERBOS
		cerr<<"Component = "<<_compNum<<" with beta = "<<beta<<" logL = "<<res<<endl;
#endif
		return -res;
	}
};


class C_evalProbMixture{
public:
  C_evalProbMixture(const tree& et,
					const sequenceContainer& sc,
					stochasticProcess* pSp,
					const int componetNumber,
					const Vdouble * weights = NULL,
					unObservableData* unObservableData_p=NULL)
		: _et(et),_sc(sc),_weights(weights),_pSp(pSp), _compNum(componetNumber)
  {
	  if(unObservableData_p)
		  _unObservableData_p = unObservableData_p->clone();
	  else
		  _unObservableData_p = NULL;
  }
  virtual ~C_evalProbMixture(){
	  if(_unObservableData_p)	delete _unObservableData_p;
  }

private:
	const tree& _et;
	const sequenceContainer& _sc;
	const Vdouble * _weights;
	stochasticProcess* _pSp;
	const int _compNum;
	unObservableData* _unObservableData_p;
public:
	MDOUBLE operator() (MDOUBLE w) {
		mixtureDistribution * pMixture = static_cast<mixtureDistribution*>(_pSp->distr());
		pMixture->setComponentWeight(w, _compNum);
		if(_unObservableData_p){
			_unObservableData_p->setLforMissingData(_et,_pSp);
		}
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,*_pSp,_weights,_unObservableData_p);
		return -res;
	}
};


/*
//the function to optimize using the conjugate Gradient algorithm
class C_evalGammaMixture {
public:
	C_evalGammaMixture(tree* pT,
					sequenceContainer* pSc,
					stochasticProcess* pSp,
					const Vdouble * weights = NULL,
					const MDOUBLE gradEps = 0.001)
		: _pTree(pT),_pSc(pSc),_pWeights(weights),_pSp(pSp), _gradEpsilon(gradEps)
	{};

		
	C_evalGammaMixture() {}

	C_evalGammaMixture& operator= (const C_evalGammaMixture &other)
	{
		_pTree = other._pTree;
		_pSc = other._pSc;
		_pWeights = other._pWeights;
		_pSp = other._pSp;
		_gradEpsilon = other._gradEpsilon;
		return *this;
	}
	
	MDOUBLE operator () (Vdouble &param){
		mixtureDistribution * pMixture = static_cast<mixtureDistribution*>(_pSp->distr());
	
		int paramNum = 0;
		for (int comp = 0; comp < pMixture->getComponentsNum(); ++comp)
		{
			pMixture->setAlpha(param[paramNum++], comp);
			pMixture->setBeta(param[paramNum++], comp);			
			pMixture->setComponentWeight(param[paramNum++], comp);			
		}
		pMixture->normalizeProbabilities();

		if (checkOutOfBounds(pMixture) == true)
			return 1000000000;
        
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(*_pTree,*_pSc,*_pSp,_pWeights);
		return -res;
	}

	void dfunc(const Vdouble &paramsIn, Vdouble& grads){
		if (paramsIn.size() != grads.size())
			errorMsg::reportError("C_evalGammaMixture::dfunc(): vectors of prameters and gradients are not the same size"); 
		Vdouble myx = paramsIn;	// temporary vector, since x is const. 

		// calc the likelihood at the current point
		MDOUBLE fa = (*this)(myx);
		
		// then calc likelihood at param+deltah for each param to approximate the derivative.
		int curParam;
		for(curParam=0; curParam < paramsIn.size(); curParam++)
		{
			myx[curParam] += _gradEpsilon;
			MDOUBLE fb = (*this)(myx);
			grads[curParam] = (fb - fa)/_gradEpsilon;
			myx[curParam] -= _gradEpsilon;
		}
	}

private:
	bool checkOutOfBounds(mixtureDistribution * pMixture) {
	for (int comp = 0; comp < pMixture->getComponentsNum(); ++comp)
	{
		if ((pMixture->getAlpha(comp) >= 15) || (pMixture->getAlpha(comp) <= 0.05))
			return true;
		if ((pMixture->getBeta(comp) >= 15) || (pMixture->getBeta(comp) <= 0.05))
			return true;
		if ((pMixture->getComponentProb(comp) > 1.0) || (pMixture->getComponentProb(comp) < 0.0))
			return true;
		}
	return false;
	}

private:
	tree* _pTree;
	sequenceContainer* _pSc;
	const Vdouble * _pWeights;
	stochasticProcess* _pSp;
	MDOUBLE _gradEpsilon; //the epsilon to calculate the gradiante
};
*/



#endif

