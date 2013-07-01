#include "optGammaMixtureLS.h"
#include "likelihoodComputation.h"
#include "numRec.h"
//#include "optimizer.h"
//#include "NRconjugateGradient.h"

#include <fstream>
#include <algorithm>
#include <ctime>
using namespace std;
using namespace likelihoodComputation;

optGammaMixtureLS::optGammaMixtureLS(stochasticProcess* pSp, const sequenceContainer& sc, const tree& inTree, MDOUBLE upperBoundAlpha/*=15.0*/, MDOUBLE upperBoundBeta/*=15.0*/,unObservableData* unObservableData_p)
{
	_pSc = &sc;
	_pTree = &inTree;
	_pSp = pSp; 
	_upperBoundAlpha = upperBoundAlpha;
	_upperBoundBeta = upperBoundBeta;
	_unObservableData_p = unObservableData_p;
}


optGammaMixtureLS::~optGammaMixtureLS()
{
}

MDOUBLE optGammaMixtureLS::optimizeParam(const int maxIterations, const MDOUBLE tol, const Vdouble * pWeights, optAlg optType)
{
	mixtureDistribution * pMixture = static_cast<mixtureDistribution*>(_pSp->distr());
	return optimizeParam(pMixture, maxIterations, tol, pWeights, optType);
}


MDOUBLE optGammaMixtureLS::optimizeParam(mixtureDistribution * pMixture, const int maxIterations, const MDOUBLE tol, const Vdouble * pWeights, optAlg optType)
{
	switch (optType) 
	{
	case ONE_DIM:
		return optimizeParamOneDim(pMixture, maxIterations, tol, pWeights);
		break;
	//case POWELL:
	//	return optimizeParamPowell(pMixture, maxIterations, tol, pWeights, pOutF);
	//	break;
	//case CONJUGATE_DERIVATIVES:
	//	return optimizeParamConjugateDeriv(pMixture, maxIterations, tol, pWeights, pOutF);
	//	break;
	default:
		errorMsg::reportError("unknown optimization algorithm in optGammaMixtureLS::optimizeParam()");
		return -1;
	}
}


//this function finds the best mixture param using a line search maximization. Each time only one parameter is optimized using the regular brent algorithm.
//CAN BE USED FOR 2 COMPONENTS ONLY (the maximization on components probabilities maximize only P1, the prob of the first component, while the prob of the second is set to 1-P1)
//total there are 5 parameters to optimize: alpha1, beta1, alpha2, beta2, and P1
MDOUBLE optGammaMixtureLS::optimizeParamOneDim(mixtureDistribution * pMixture, const int maxIterations, const MDOUBLE tol, const Vdouble * pWeights)
{
	MDOUBLE lowerBound = 0.0;	
	
	MDOUBLE newL = VERYSMALL; //newL is the LL after a single param optimization.
	//MDOUBLE curL = VERYSMALL; //the current LL.
	MDOUBLE curL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(*_pTree,*_pSc,*_pSp,pWeights,_unObservableData_p); //the current LL.
	MDOUBLE prevIterL = VERYSMALL; //The LL of the previous iteration. the loop quit if the increase in LL between iterations is smaller than tol
	MDOUBLE bestA=0, bestB=0, bestW = 0;
			
	for (int it = 0; it < maxIterations; ++it)
	{
		//prevIterL = newL;
		prevIterL = curL;

		for (int comp = 0; comp < pMixture->getComponentsNum(); ++comp)
		{
			//optimize alpha 
			MDOUBLE oldAlpha = pMixture->getAlpha(comp);
			newL = -brent(lowerBound,oldAlpha, _upperBoundAlpha, C_evalAlphaMixture(*_pTree,*_pSc,_pSp,comp,pWeights,_unObservableData_p), tol, &bestA);
			if (newL < curL) 
			{
				//the Likelihood wend down
				pMixture->setAlpha(oldAlpha, comp);
				if(_unObservableData_p){
					_unObservableData_p->setLforMissingData(*_pTree,_pSp);				
				}
				LOG(5, <<"likelihood went down in optGammaMixtureLS::optimizeParam()"<<endl<<"old L= "<<curL<<" newL = "<<newL<<endl);
			}
			else
			{
				pMixture->setAlpha(bestA, comp);
				if(_unObservableData_p){
					_unObservableData_p->setLforMissingData(*_pTree,_pSp);				
				}
				curL = newL;
				LOG(7, <<"iteration: "<<it<<" Optimize alpha comp"<<comp<<" new Likelihood = "<<curL<<endl);
			}
		
			//optimize beta 
			MDOUBLE oldBeta = pMixture->getBeta(comp);
			newL = -brent(lowerBound,oldBeta,_upperBoundBeta, C_evalBetaMixture(*_pTree,*_pSc,_pSp,comp,pWeights,_unObservableData_p), tol, &bestB);
			if (newL < curL)
			{
				//the Likelihood wend down
				pMixture->setBeta(oldBeta, comp);
				if(_unObservableData_p){
					_unObservableData_p->setLforMissingData(*_pTree,_pSp);				
				}
				LOG(5, <<"likelihood went down in optGammaMixtureLS::optimizeParam()"<<endl<<"old L= "<<curL<<" newL = "<<newL<<endl);
			}
			else
			{
				pMixture->setBeta(bestB, comp);
				if(_unObservableData_p){
					_unObservableData_p->setLforMissingData(*_pTree,_pSp);				
				}
				curL = newL;
				LOG(7, <<"iteration: "<<it<<" Optimize beta comp"<<comp<<" new Likelihood = "<<curL<<endl);
			}
			//optimize components probability. 
            if (pMixture->getComponentsNum() == 1)
				continue;
		
			MDOUBLE upperBound = 0.0;
			MDOUBLE lowerBound = 1.0;
			MDOUBLE oldWeight = pMixture->getComponentWeight(comp);
			newL = -brent(lowerBound, oldWeight, upperBound, C_evalProbMixture(*_pTree,*_pSc, _pSp, comp, pWeights), tol, &bestW);
			if (newL < curL)
			{
				//the Likelihood wend down
				pMixture->setComponentWeight(oldWeight, comp);
				if(_unObservableData_p){
					_unObservableData_p->setLforMissingData(*_pTree,_pSp);				
				}
				LOG(5, <<"likelihood went down in optGammaMixtureLS::optimizeParam()"<<endl<<"old L= "<<curL<<" newL = "<<newL<<endl);
			}
			else
			{
				pMixture->setComponentWeight(bestW, comp);
				if(_unObservableData_p){
					_unObservableData_p->setLforMissingData(*_pTree,_pSp);				
				}
				curL = newL;
				LOG(7, <<"iteration: "<<it<<" Optimize Prob"<<" new Likelihood = "<<curL<<endl);
			}
		}
		pMixture->normalizeProbabilities();	// why again ???
		printIter(pMixture, it, curL);
		if (curL < prevIterL + tol){
			//if(_unObservableData_p){
			//	_unObservableData_p->setLforMissingData(*_pTree,_pSp);				
			//}
			return max(curL,prevIterL); // not to reduce likelihood
		}
	}
	return curL;
}



/*
//this function uses a line search maximization. The difference is that it does not use the naive method (optimize each parameter seperatly untill convergence)
//but uses Powel's quadratically convergent method (Numerical Recipes pp 420).
//CAN BE USED FOR 2 COMPONENTS ONLY (the maximization on components probabilities maximize only P1, the prob of the first component, while the prob of the second is set to 1-P1)
//total there are 5 parameters to optimize: alpha1, beta1, alpha2, beta2, and P1
MDOUBLE optGammaMixtureLS::optimizeParamPowell(mixtureDistribution* pMixture, const int maxIterations, const MDOUBLE tol, const Vdouble * pWeights, ofstream* pOutF)
{
	if (pMixture->getComponentsNum() == 1)
		return optimizeParam1CompPowel(pMixture, maxIterations, tol, pWeights, pOutF);
	else return optimizeParamManyCompPowel(pMixture, maxIterations, tol, pWeights, pOutF);
}


MDOUBLE optGammaMixtureLS::optimizeParam1CompPowel(mixtureDistribution * pMixture, const int maxIterations, const MDOUBLE tol, const Vdouble * pWeights, ofstream* pOutF)
{
	tree tree1(*_pTree);
	sequenceContainer sc1(*_pSc);

	C_evalGammaMixture optPowell(&tree1, &sc1, _pSp, NULL);
	optimizer<C_evalGammaMixture> opt(optPowell);
	Vdouble param(2);
	param[0] = pMixture->getAlpha(0);
	param[1] = pMixture->getBeta(0);

	MDOUBLE res = opt.findmin(param);
	return res;
}

MDOUBLE optGammaMixtureLS::optimizeParamManyCompPowel(mixtureDistribution * pMixture, const int maxIterations, const MDOUBLE tol, const Vdouble * pWeights, ofstream* pOutF)
{
	tree tree1(*_pTree);
	sequenceContainer sc1(*_pSc);

    Vdouble param(pMixture->getComponentsNum() * 3 - 1);
	int paramNum = 0;
	for (int comp = 0; comp < pMixture->getComponentsNum(); ++comp)
	{
        param[paramNum++] = pMixture->getAlpha(comp);
        param[paramNum++] = pMixture->getBeta(comp);
		param[paramNum++] = pMixture->getComponentWeight(comp);
	}
	C_evalGammaMixture optPowell(&tree1, &sc1, _pSp, NULL);
	optimizer<C_evalGammaMixture> opt(optPowell);
	MDOUBLE res = opt.findmin(param);
	cerr <<"optimized Powell result = "<< res<<endl;
	return res;
}
*/

/*
MDOUBLE optGammaMixtureLS::optimizeParamConjugateDeriv(
	mixtureDistribution * pMixture, const int maxIterations, 
	const MDOUBLE tol, const Vdouble * pWeights, ofstream* pOutF)
{
	tree tree1(*_pTree);
	sequenceContainer sc1(*_pSc);

    Vdouble param(pMixture->getComponentsNum() * 3);
	int paramNum = 0;
	int comp;
	for (comp = 0; comp < pMixture->getComponentsNum(); ++comp)
	{
        param[paramNum++] = pMixture->getAlpha(comp);
        param[paramNum++] = pMixture->getBeta(comp);
		param[paramNum++] = pMixture->getComponentWeight(comp);
	}
	C_evalGammaMixture func(&tree1, &sc1, _pSp, pWeights);
	NRconjugateGradient<C_evalGammaMixture> opt;
	if (pOutF != NULL)
	{
		*pOutF <<endl<<endl<<"starting NRconjugateGradient optimization..."<<endl;
		printIter(pMixture, 0, 0.0, pOutF);
	}

	MDOUBLE res = opt.findmin(param, &func, tol);

	paramNum = 0;
	for (comp = 0; comp < pMixture->getComponentsNum(); ++comp)
	{
		pMixture->setAlpha(param[paramNum++], comp);
		pMixture->setBeta(param[paramNum++], comp);			
		pMixture->setComponentWeight(param[paramNum++], comp);			
	}
	pMixture->normalizeProbabilities();
	if (pOutF != NULL)
	{
		*pOutF <<endl<<endl<<"after NRconjugateGradient optimization"<<endl;
		printIter(pMixture, 0, res, pOutF);
	}
	cerr <<"optimized Conjugate Deriv result = "<< res<<endl;
	return res;
}
*/


void optGammaMixtureLS::printIter(const mixtureDistribution * pMixture, const int it, const MDOUBLE curL)  
{
	LOG(4,<< "iter " << it <<": cur likelihood= " << curL <<endl);
	for (int k = 0; k < pMixture->getComponentsNum(); ++k)
	{
		LOG(4, << "comp="<<k<<" Alp/Beta= "<<pMixture->getAlpha(k)/pMixture->getBeta(k)<<" alpha= "<<pMixture->getAlpha(k) << " beta= " <<pMixture->getBeta(k)<<" Prob= "<<pMixture->getComponentProb(k)<<endl);  
	}
}
