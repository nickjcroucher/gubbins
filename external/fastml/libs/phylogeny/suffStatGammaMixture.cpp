#include "suffStatGammaMixture.h"
#include "mixtureDistribution.h"
#include "computePijComponent.h"
#include "likelihoodComputation.h"
#include "gammaUtilities.h"
#include "uniDistribution.h"


#include <cmath>
#include <fstream>
using namespace likelihoodComputation;


suffStatGammaMixture::suffStatGammaMixture(const stochasticProcess& cur_sp, const sequenceContainer& sc, const tree& inTree)
{
	_pSp = &cur_sp;
	_pSc = &sc;
	_pTree = &inTree;
}

suffStatGammaMixture::~suffStatGammaMixture()
{
}


void suffStatGammaMixture::allocatePlaceForSuffStat() {
	mixtureDistribution* pMixture = static_cast<mixtureDistribution*>(_pSp->distr());
	int componentNum = pMixture->getComponentsNum();
	_MkVec.clear();
	_MkVec.resize(componentNum, 0);
	_AkVec.clear();
	_AkVec.resize(componentNum, 0);
	_BkVec.clear();
	_BkVec.resize(componentNum, 0);
}

void suffStatGammaMixture::computePijForEachComponent(vector<computePijGam>& cpgVec,
													  vector<stochasticProcess>& spVec) {
	mixtureDistribution* pMixture = static_cast<mixtureDistribution*>(_pSp->distr());
	int componentNum = pMixture->getComponentsNum();
	for (int comp = 0; comp < componentNum; ++comp) {
		//create a local sp so to compute likelihoods of this component only
		stochasticProcess compSp(pMixture->getComponent(comp), _pSp->getPijAccelerator());
		cpgVec[comp].fillPij(*_pTree, compSp);
		spVec.push_back(compSp);
	}
}

void suffStatGammaMixture::computeStatistics()
{
	///////////////as in getTreeLikelihoodAllPosAlphTheSame
	//computePijGam pi;
	//pi.fillPij(*_pTree, *_pSp);
	//MDOUBLE res =0;
	//doubleRep LofPos;
	//int k;
	//for (k=0; k < _pSc->seqLen(); ++k) 
	//{
	//	doubleRep tmp=0;
	//	for (int i=0; i < _pSp->categories();++i) 
	//	{
	//		tmp += getLofPos(k, *_pTree, *_pSc, pi[i], *_pSp) *  _pSp->ratesProb(i);
	//	}
	//	LofPos = tmp;
	//	res += log(LofPos);
	//}
	//////////////////////////////////////////////

	//mixtureDistribution* pMixture = static_cast<mixtureDistribution*>(_pSp->distr());
	//int componentNum = pMixture->getComponentsNum();
	//MDOUBLE res2 = 0.0;
	//vector<computePijGam> cpgVec(componentNum);
	//vector<stochasticProcess> spVec;
	//
	//for (int comp = 0; comp < componentNum; ++comp) {
	//	//create a local sp so to compute likelihoods of this component only
	//	stochasticProcess compSp(pMixture->getComponent(comp), _pSp->getPijAccelerator());
	//	cpgVec[comp].fillPij(*_pTree, compSp);
	//	spVec.push_back(compSp);
	//}
	//
	//for (int pos = 0; pos < _pSc->seqLen(); ++pos)
	//{
	//	int comp;
	//	for (comp = 0; comp < componentNum; ++comp)
	//	{
	//		const generalGammaDistribution* pDist = pMixture->getComponent(comp);
	//		for (int cat=0; cat < pDist->categories(); ++cat) 
	//		{
	//			MDOUBLE tmp = pDist->ratesProb(cat) * getLofPos(pos, *_pTree, *_pSc, cpgVec[comp][cat], *_pSp); 
	//			res2 += log(tmp);
	//		}
	//	}
	//}
	//////////////////////////////////////////////
	allocatePlaceForSuffStat();
	mixtureDistribution* pMixture = static_cast<mixtureDistribution*>(_pSp->distr());
	int componentNum = pMixture->getComponentsNum();

	//compute Pij for each component
	vector<computePijGam> cpgVec(componentNum);
	vector<stochasticProcess> spVec;
	computePijForEachComponent(cpgVec,spVec);


	//compute statistics: M_k, A_k, B_k
	//Here we sum over all positions.
	//go over all positions [pos] and compute for each component [k]: M_k(pos), E[R]_k(pos), E[logR]_k(pos)
	//Then compute A_k and B_k for that position.
	for (int pos = 0; pos < _pSc->seqLen(); ++pos)
	{
		MDOUBLE sumAllComponents = 0.0;
		Vdouble MkPosVec(componentNum, 0.0); //the contribution of position pos to the M_K statistic
		Vdouble Exp_RkVec(componentNum, 0.0);
		Vdouble Exp_LogRkVec(componentNum, 0.0);
		int comp;
		for (comp = 0; comp < componentNum; ++comp)
		{
			// here we compute P(H[i]=k, data| cur_mixtureDistribution)
			//P(H[i]=k, data| teta) = P(H[i]=k)* (sum_over_all_categories{P(data|r)P(r))
			///////////////////////////
			const generalGammaDistribution* pDist = pMixture->getComponent(comp);
			MDOUBLE Exp_Rk, Exp_LogRk, sum;
			Exp_Rk = Exp_LogRk = sum = 0.0;
			for (int cat=0; cat < pDist->categories(); ++cat) 
			{
				MDOUBLE LofP = convert(likelihoodComputation::getLofPos(pos, *_pTree, *_pSc, cpgVec[comp][cat], spVec[comp]));
				MDOUBLE Pr = pDist->ratesProb(cat) * LofP; 
				sum += Pr;
				Exp_RkVec[comp] += Pr * pDist->rates(cat);
				Exp_LogRkVec[comp] += Pr * log(pDist->rates(cat));
			}
			MkPosVec[comp] = sum;
			sumAllComponents += MkPosVec[comp] * pMixture->getComponentProb(comp);;
		}

		for (comp = 0; comp < componentNum; ++comp)
		{
			MDOUBLE factor = pMixture->getComponentProb(comp)/ sumAllComponents;
			_MkVec[comp] += factor* MkPosVec[comp] ;
			_AkVec[comp] += factor * Exp_RkVec[comp];
			_BkVec[comp] += factor * Exp_LogRkVec[comp];
		}
	}// end of loop over positions
	spVec.clear();
	cpgVec.clear();
}


#include "uniformDistribution.h"
void suffStatGammaMixture::plotStatistics(ofstream& outFile)
{
	mixtureDistribution* pMixture = static_cast<mixtureDistribution*>(_pSp->distr());
	if (pMixture->getComponentsNum() != 1)
		errorMsg::reportError("Sorry, I plot only 1 component");

	outFile <<"R"<<"\t"<<"Postr"<<"\t"<<"Er"<<"\t"<<"Elog_r"<<endl;
	const generalGammaDistribution* pDist = pMixture->getComponent(0);
	int numCat = 200, maxR = 10;
	uniformDistribution uniDist(numCat, 0, maxR);
	/////////calc the prior of each interval
	Vdouble priorProbs(uniDist.categories());
	MDOUBLE upperP, lowerP = 0;  
	for (int i = 0; i<uniDist.categories();++i)
	{
		upperP = pDist->getCumulativeProb(uniDist.getBorder(i+1));
		priorProbs[i] = upperP - lowerP;
		lowerP = upperP;
	}

	distribution * pUni =  new uniDistribution;

	stochasticProcess uniSp(pUni, _pSp->getPijAccelerator()); 
	//loop over all r
	for (int ri=0; ri < uniDist.categories(); ++ri) 
	{
		MDOUBLE Exp_R = 0.0;
		MDOUBLE Exp_LogR = 0.0;
		MDOUBLE PosteriorR = 0.0;
		MDOUBLE rate = uniDist.rates(ri);
		if (rate == 0.0)
			rate = 0.000001;
		
		//Here we sum over all positions.
		//go over all positions [pos] and compute: PosrteriorR(=P(D|r)*P(r)), E[R]_k(pos), E[logR]_k(pos)
		for (int pos = 0; pos < _pSc->seqLen(); ++pos)
		{
			MDOUBLE PrPos = priorProbs[ri] * convert(likelihoodComputation::getLofPos(pos, *_pTree, *_pSc, uniSp, rate)); 
			PosteriorR += PrPos;
			Exp_R += PrPos * rate;
			Exp_LogR += PrPos * log(rate);

		}

		outFile <<rate<<"\t"<<PosteriorR<<"\t"<<Exp_R<<"\t"<<Exp_LogR<<endl;
	}

	delete pUni;
}


MDOUBLE suffStatGammaMixture::computeQ2() 
{
	MDOUBLE res=0;

	return res;
}



MDOUBLE suffStatGammaMixture::computeQ() 
{
	mixtureDistribution* pMixture = static_cast<mixtureDistribution*>(_pSp->distr());
	MDOUBLE res = 0.0;
	//////////////////////////////////
	MDOUBLE res2 = 0.0;
	int compNum = pMixture->getComponentsNum();
	///////////////////////////////////
	for (int comp = 0;comp < compNum ; ++comp)
	{
		MDOUBLE P_k = pMixture->getComponentProb(comp);
		MDOUBLE alpha_k = pMixture->getAlpha(comp);
		MDOUBLE beta_k = pMixture->getBeta(comp);
		MDOUBLE first = _MkVec[comp] * log(P_k);
		MDOUBLE second = _MkVec[comp] * alpha_k*log(beta_k);
		MDOUBLE third = -_MkVec[comp] * gammln(alpha_k);
		MDOUBLE fourth = -_AkVec[comp]*beta_k;
		MDOUBLE fifth = _BkVec[comp]*(alpha_k-1.0);
		res += _MkVec[comp] * (log(P_k) + alpha_k*log(beta_k) - gammln(alpha_k)) 
			   - (_AkVec[comp]*beta_k) 
			   + _BkVec[comp]*(alpha_k-1);
		////////////////////////////////////
	}
	res2 = computeQ2();
	return res;
}
