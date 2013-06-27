#include "GamMixtureOptimizer.h"
#include "someUtil.h"
#include "optGammaMixtureEM.h"
#include "optGammaMixtureLS.h"

#include <fstream>
#include <algorithm>
#include <ctime>
using namespace std;



GamMixtureOptimizer::GamMixtureOptimizer(stochasticProcess* pSp, const sequenceContainer& sc, const tree& inTree, unObservableData* unObservableData_p)
{
	_pSc = &sc;
	_pTree = &inTree;
	_pSp = pSp;
	_unObservableData_p = unObservableData_p;
	_tolOptSpecific = 0.001;

}


GamMixtureOptimizer::~GamMixtureOptimizer()
{
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//findBestParamManyStarts: Finds the best gammaMixture from many starting points.
//The function starts form few starting points. 
//For each point it tries to optimize the likellihood doing only a small number of iterations.
//It then picks the best points (highest likelihood) and continue the maximization for these points only.
//This can be repeated a number of times, each cycle with a different optimization algorithm.
//The best gammaMixture is stored in _sp and the best likelihood is returned.
//input Parameters:
//pointsNum: a vector with the number of points to peformed the current cycle of optimization.
//iterNum: the number of iterations to perform in each cycle.
//OptAlgs: the optimization algorithm to be performed in each cycle.
//tol        = for determining convergence in the maximization process. 
MDOUBLE GamMixtureOptimizer::findBestParamManyStarts(const Vint pointsNum, const Vint iterNum, const vector<OptimAlg> OptAlgs, const Vdouble tols, const Vdouble * pWeights, ofstream* pOutF/*= NULL*/)
{
	//make sure that the number of points in each cycle is not bigger than the previous cycle.
	int i;	
	for (i = 0; i < pointsNum.size()-1; ++i)
	{
		if (pointsNum[i] < pointsNum[i+1])
			errorMsg::reportError("input error in GamMixtureOptimizer::findBestParamManyStarts()");
	}

	//create starting distributions 
	vector<mixtureDistribution*> distVec;
	const mixtureDistribution * pMixture = getMixtureDist();
	for (i = 0; i < pointsNum[0]; ++i)
	{
		//the first distribution will be the current one
		if (i == 0)
			distVec.push_back(new mixtureDistribution(*pMixture)); 
		else
			distVec.push_back(new mixtureDistribution(pMixture->getComponentsNum(), pMixture->categoriesForOneComponent(), LAGUERRE, 15, 15)); 
	}

	//make a small number of iterations for all random starts 
	int numOfOptCycles = pointsNum.size();
	Vdouble likelihoodVec;
	for (i = 0; i < numOfOptCycles; ++i)
	{
		if (i != 0)
		{
			vector<mixtureDistribution*> tmpDistVec(0);
			//sort results and continue optimization only with the best (pointsNum[i]) points
			Vdouble sortedL = likelihoodVec;
			sort(sortedL.begin(),sortedL.end());
			MDOUBLE threshold = sortedL[sortedL.size()- pointsNum[i]];
			for (int j = 0; j < likelihoodVec.size(); ++j)
			{
				if (likelihoodVec[j] >= threshold) 
					tmpDistVec.push_back(distVec[j]);
				else
					delete distVec[j];
			}
			distVec.clear();
			distVec = tmpDistVec;
		} 
		
		likelihoodVec.clear();
		likelihoodVec.resize(pointsNum[i]); 
		int c; 		
		for (c = 0; c < pointsNum[i]; ++c)
		{
			cerr <<"optimizing point " <<c<<endl;
			MDOUBLE ll = optimizeParam(distVec[c], iterNum[i], OptAlgs[i], tols[i], pWeights, pOutF);
			cerr<<"pointi: "<<c<<"  likelihood = "<<ll<<endl;
			likelihoodVec[c] = ll;
		}
	}

	Vdouble sortedL = likelihoodVec;
	sort(sortedL.begin(),sortedL.end());
	MDOUBLE bestL = sortedL[likelihoodVec.size() - 1];
	for (i = 0; i < likelihoodVec.size(); ++i)
	{
		if (bestL == likelihoodVec[i]) 
		{
			_pSp->setDistribution(distVec[i]);
		}
		delete distVec[i];
	}	
	distVec.clear();
	return bestL;
}

MDOUBLE GamMixtureOptimizer::findBestParam(const OptimAlg alg, const int maxIterations, const MDOUBLE tol, const Vdouble * pWeights, ofstream* pOutF)
{
	mixtureDistribution* pInDistribution = static_cast<mixtureDistribution*>(_pSp->distr());
	return optimizeParam(pInDistribution, maxIterations, alg, tol, pWeights, pOutF);
}


MDOUBLE GamMixtureOptimizer::optimizeParam(mixtureDistribution* pInDistribution, const int maxIterations, const OptimAlg alg, const MDOUBLE tol, const Vdouble * pWeights, ofstream* pOutF)
{
	MDOUBLE res = 0.0;
	switch (alg)
	{
	case EM: {
		optGammaMixtureEM emOpt(*_pSp, *_pSc, *_pTree); 
		res = emOpt.optimizeParam(pInDistribution, maxIterations, tol, _tolOptSpecific, pOutF);
		break;
		}
	case ONE_DIM: {
		optGammaMixtureLS lsOpt(_pSp, *_pSc, *_pTree,MAXIMUM_ALPHA_PARAM,MAXIMUM_BETA_PARAM,_unObservableData_p);
		res = lsOpt.optimizeParam(pInDistribution, maxIterations, tol, pWeights, optGammaMixtureLS::ONE_DIM);
		MDOUBLE resRecompute = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(*_pTree,*_pSc,*_pSp,pWeights,_unObservableData_p);
		if(!DEQUAL(res,resRecompute)){
			LOGnOUT(3,<<"--- error: different likelihood after GamMixtureOptimizer::optimizeParam,diff= "<<res-resRecompute <<"\n");
		}
		break;
		}
	//case TX_CONJUGATE_DERIVATIVES:
	//	{
	//	txGamMixtureOptimizer txOpt(_pSp, *_pSc, *_pTree);
	//	txOpt.setOptimizationParameters(tol, _tolOptSpecific, _tolOptSpecific, _tolOptSpecific);
	//	res = txOpt.optimizeParam(pInDistribution, maxIterations, pWeights, alg, pOutF);
	//	break;
	//	}
	//case NR_CONJUGATE_DERIVATIVES:
	//	{
	//	optGammaMixtureLS opt(_pSp, *_pSc, *_pTree);
	//	res = opt.optimizeParam(pInDistribution, maxIterations, tol, pWeights, optGammaMixtureLS::CONJUGATE_DERIVATIVES, pOutF);
	//	break;
	//	}
	default:
		errorMsg::reportError("unknown optimization algorithm in GamMixtureOptimizer::optimizeParam()");
	}
	return res;
}
