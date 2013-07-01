#include "optGammaMixtureEM.h"
#include "likelihoodComputation.h"
#include "numRec.h"
#include "uniDistribution.h"

#include <fstream>
#include <algorithm>
#include <ctime>
using namespace std;
using namespace likelihoodComputation;

optGammaMixtureEM::optGammaMixtureEM(const stochasticProcess& cur_sp, const sequenceContainer& sc, const tree& inTree)
{
	_pSc = &sc;
	_pTree = &inTree;
	_pSp = new stochasticProcess(cur_sp); 
}

optGammaMixtureEM::~optGammaMixtureEM()
{
	if (_pSp != NULL)
	{
		delete _pSp;
		_pSp = NULL;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//findBestParamManyStarts: Finds the best gammaMixture from many starting points.
//The function starts form few starting points. 
//For each point it tries to optimize the likellihood doing only a small number of iterations.
//It then picks the best points (highest likelihood) and continue the maximization for these points only.
//The best gammaMixture is stored in _sp and the best likelihood is returned.
//input Parameters:
//startPointsNum = the number of starting points.
//bestStartsNum  = the number of best points to continue with the full optimization.
//startIter      = the number of iterations to perform with all starting points.
//maxIterations  = the maximum number of iterations to continue with the best points
//epsilon        = for determining convergence in the maximization process. 
MDOUBLE optGammaMixtureEM::findBestParamManyStarts(const int startPointsNum, const int bestStartsNum, const int startIter, const int maxIterations, const MDOUBLE epsilon, const MDOUBLE epsilomQopt, ofstream* pOutF)
{
	vector<mixtureDistribution> distVec;
	Vdouble likelihoodVec(startPointsNum);
	mixtureDistribution * pMixture = static_cast<mixtureDistribution*>(_pSp->distr());
	//create starting distributions
	int i;
	for (i = 0; i < startPointsNum; ++i)
	{
		//the first distribution will be the current one
		if (i == 0)
			distVec.push_back(*pMixture); 
		else
			distVec.push_back(mixtureDistribution(pMixture->getComponentsNum(), pMixture->categoriesForOneComponent(), LAGUERRE, 15, 15)); 
	}

	//make a small number of iterations for all random starts 
	for (i = 0; i < distVec.size(); ++i)
	{
		likelihoodVec[i] = optimizeParam(&distVec[i], startIter, epsilon, epsilomQopt, pOutF);
	}

	//sort results and make full optimization only on the best starts
	Vdouble sortedL = likelihoodVec;
	sort(sortedL.begin(),sortedL.end());
	MDOUBLE threshold = sortedL[sortedL.size()- bestStartsNum];
	MDOUBLE bestL = sortedL[0];
	int bestDistNum = 0;
	for (i = 0; i < distVec.size(); ++i)
	{
		if (likelihoodVec[i] >= threshold) 
		{
			MDOUBLE newL = optimizeParam(&distVec[i], maxIterations, epsilon, epsilomQopt, pOutF);
			if (newL > bestL)
			{
				bestL = newL;
				bestDistNum = i;
			}
		}
	}
	_pSp->setDistribution(&distVec[bestDistNum]);
	distVec.clear();
	return bestL;
}


MDOUBLE optGammaMixtureEM::optimizeParam(mixtureDistribution* pInDistribution, const int maxIterations, const MDOUBLE epsilon, const MDOUBLE epsilomQopt, ofstream* pOutF)
{
	stochasticProcess inSp(pInDistribution, _pSp->getPijAccelerator());
	MDOUBLE curL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(*_pTree, *_pSc, inSp, NULL);

	/////compute piHomPos as in getTreeLikelihoodAllPosAlphTheSame
	//computePijGam pi;
	//pi.fillPij(*_pTree, inSp);
	//MDOUBLE res =0;
	//doubleRep LofPos;
	//int k;
	//for (k=0; k < _pSc->seqLen(); ++k) 
	//{
	//	doubleRep tmp=0;
	//	for (int i=0; i < inSp.categories();++i) 
	//	{
	//		tmp += getLofPos(k, *_pTree, *_pSc, pi[i], inSp)* inSp.ratesProb(i);
	//		/*MDOUBLE Pr = pDist->ratesProb(cat) * likelihoodComputation::getLofPos(pos, *_pTree, *_pSc, cpgVec[comp][cat], spVec[comp]); */
	//	}
	//	LofPos = tmp;
	//	res += log(LofPos);
	//}
	//




	
	
	//int componentNum = pInDistribution->getComponentsNum();
	////compute Pij for each component
	//vector<computePijGam> cpgVec(componentNum);
	//vector<stochasticProcess> spVec;
	//for (int comp = 0; comp < componentNum; ++comp) {
	//	//create a local sp so to compute likelihoods of this component only
	//	stochasticProcess compSp(pInDistribution->getComponent(comp), _pSp->getPijAccelerator());
	//	cpgVec[comp].fillPij(*_pTree, compSp);
	//	spVec.push_back(compSp);
	//}



	//for (int pos = 0; pos < _pSc->seqLen(); ++pos)
	//{
	//	int comp;
	//	for (comp = 0; comp < componentNum; ++comp)
	//	{
	//		const generalGammaDistribution* pDist = pInDistribution->getComponent(comp);
	//		for (int cat=0; cat < pDist->categories(); ++cat) 
	//		{
	//			doubleRep LofPos = likelihoodComputation::getLofPos(pos, *_pTree, *_pSc, cpgVec[comp][cat], spVec[comp]); 
	//			L2 += log(LofPos);
	//		}
	//	}
	//}


	
	if (maxIterations == 0)
	{
		return curL;
		LOG(4,<<endl<<endl<<"starting Gamma Mixture EM optimization..."<<endl);
		printIter(inSp, 0, curL);
	}

	MDOUBLE newL = curL;
	int it;
	for (it = 0; it < maxIterations; ++it)
	{
		stochasticProcess oldSp(inSp);
		maximizeGammaParam(&inSp, epsilomQopt);
		newL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(*_pTree, *_pSc, inSp, NULL);
		if (newL < curL + epsilon)
		{
			//the improvemnt in Likelihood is smaller than epsilon
			if (newL < curL)
			{	//ERROR - L went Down!
				cerr<<"likelihood went down!"<<endl<<"oldL = "<<curL<<" newL= "<<newL<<" Diff= "<<newL-curL<<endl;
				if (pOutF != NULL) *pOutF <<"likelihood went down!"<<endl<<"oldL = "<<curL<<" newL= "<<newL<<endl;
				*pInDistribution = *(static_cast<mixtureDistribution*>(oldSp.distr()));
				if (pOutF != NULL) *pOutF <<"after Gamma Mixture EM optimization..."<<endl;
				printIter(inSp, it, curL);
				return curL;
			}
			else
			{
				cerr<<"converged!"<<endl; 
				*pInDistribution = *(static_cast<mixtureDistribution*>(inSp.distr()));
				if (pOutF != NULL) *pOutF <<"after Gamma Mixture EM optimization..."<<endl;
				printIter(inSp, it, newL);
				return newL;
			}
		}
		cerr << "iter " << it <<": cur likelihood= " << curL <<" new likelihood= " << newL <<endl;
		curL = newL;
	}
	
	*pInDistribution = *(static_cast<mixtureDistribution*>(inSp.distr()));
	if (pOutF != NULL) *pOutF <<"after Gamma Mixture EM optimization..."<<endl;
	printIter(inSp, it, newL);
	return newL;
}


void optGammaMixtureEM::maximizeGammaParam(stochasticProcess* pNewSp, MDOUBLE accuracyRtbis)
{
	suffStatGammaMixture stats(*pNewSp, *_pSc, *_pTree);
	stats.computeStatistics();
	//cerr << "Q BEFORE IS: " << stats.computeQ()<<endl;
	maximizeGammaParam(stats, pNewSp, accuracyRtbis);
	//cerr << "Q AFTER IS: " << stats.computeQ()<<endl;
}

void optGammaMixtureEM::maximizeGammaParam(const suffStatGammaMixture & stats,
										   stochasticProcess* pNewSp, MDOUBLE accuracyRtbis)
{
	MDOUBLE upperBoundAlpha = 15.0;
	mixtureDistribution * pMixture = static_cast<mixtureDistribution*>(pNewSp->distr());
	int numComponents = pMixture->getComponentsNum();
	Vdouble compProb(numComponents), alphaVec(numComponents), betaVec(numComponents);
	for (int k = 0; k < numComponents; ++k)
	{
		alphaVec[k] = findBestAlpha(stats, k, accuracyRtbis, upperBoundAlpha);
		betaVec[k] = alphaVec[k] * (stats.getMk(k) / stats.getAk(k)); 
		compProb[k] = stats.getMk(k) / _pSc->seqLen();
	}
	pMixture->setMixtureParameters(alphaVec, betaVec, compProb);
}

void optGammaMixtureEM::printIter(const stochasticProcess& inSp, const int it, const MDOUBLE curL)  
{
	LOG(4, << "iter " << it <<": cur likelihood= " << curL <<endl);
	mixtureDistribution * pMixture = static_cast<mixtureDistribution*>(inSp.distr());
	for (int k = 0; k < pMixture->getComponentsNum(); ++k)
	{
		LOG(4, << "comp="<<k<<" Alp/Beta= "<<pMixture->getAlpha(k)/pMixture->getBeta(k)<<" alpha= "<<pMixture->getAlpha(k) << " beta= " <<pMixture->getBeta(k)<<" Prob= "<<pMixture->getComponentProb(k)<<endl);  
	}
}


//findBestAlpha: this function finds the alpha which is the root of the function C_evalAlphaEM().
//BUT - if there is no root in the range (lowerBoundAlpha, upperBoundAlpha) 
//or - the root is higher than upperBoundAlpha - the function returns upperBoundAlpha
MDOUBLE optGammaMixtureEM::findBestAlpha(const suffStatGammaMixture& stats, const int compNum, const MDOUBLE accuracyRtbis, const MDOUBLE upperBoundAlpha) const 
{
	MDOUBLE res = upperBoundAlpha;
	MDOUBLE lowerBoundAlpha = MINIMUM_ALPHA_PARAM;
	MDOUBLE upperRange = upperBoundAlpha;
	MDOUBLE lowerRange = lowerBoundAlpha;	
	bool haveRoot = zbrac(C_evalAlphaEM(stats, compNum), lowerRange, upperRange);
	if (haveRoot == true)
		res = rtbis(C_evalAlphaEM(stats, compNum), lowerRange, upperRange, accuracyRtbis); ;
	if (res > upperBoundAlpha)
		res = upperBoundAlpha;
	else if (res < lowerBoundAlpha)
		res = lowerBoundAlpha;
	return res;
}


void optGammaMixtureEM::checkEntropy(stochasticProcess & oldSp, stochasticProcess & newSp)
{
	//the entropy is 
	//sigma_r P(r|D,oldSp)*log(P(r|D,oldSp) / P(r|D,newSp))
	//VVdouble posteriorBefore,posteriorAfter ; 
	//likelihoodComputation::getPosteriorOfRates(*_pTree, *_pSc, oldSp, posteriorBefore, NULL);
	//likelihoodComputation::getPosteriorOfRates(*_pTree, *_pSc, newSp, posteriorAfter, NULL);


	//MDOUBLE entropyAll = 0.0;
	//MDOUBLE secondTerm= 0.0;
	//for (int pos = 0; pos < _pSc->seqLen(); ++pos)
	//{
	//	MDOUBLE entropyPos = 0.0;	
	//	for (int cat = 0; cat < oldSp.categories(); ++cat)
	//	{
	//		entropyPos += posteriorBefore[pos][cat] * log(posteriorBefore[pos][cat] / posteriorAfter[pos][cat]);
	//		secondTerm += posteriorBefore[pos][cat] * log(posteriorAfter[pos][cat]);
	//	}
	//	entropyAll += entropyPos;
	//	//cerr <<"Pos Entropy = "<<entropyPos<<endl;
	//}
	//cerr <<endl<<endl<<endl;
	//cerr <<"All Entropy = "<<entropyAll<<endl;


	//calculating Q
	//MDOUBLE QAll = 0.0;
	//for (int pos = 0; pos < _pSc->seqLen(); ++pos)
	//{
	//	MDOUBLE QPos = 0.0;	
	//	for (int cat = 0; cat < oldSp.categories(); ++cat)
	//	{
	//		stochasticProcess localSp(&uniDistribution(), newSp.getPijAccelerator()); 
	//		MDOUBLE rate = newSp.rates(cat);
	//		MDOUBLE L_after = likelihoodComputation::getLofPos(pos, *_pTree, *_pSc, localSp, rate);
	//		QPos += posteriorBefore[pos][cat] * log(L_after * newSp.ratesProb(cat));
	//	}
	//	QAll += QPos;
	//	//cerr <<"Pos Q = "<<QPos<<endl;
	//}
	//cerr <<endl<<endl<<endl;
	//cerr <<"Q ALL= "<<QAll<<endl;
	//cerr <<"secondTerm = "<<secondTerm<<endl;
	
}
