#include "mixtureDistribution.h"
#include "generalGammaDistributionLaguerre.h"
#include "talRandom.h"
#include "someUtil.h"
#include "errorMsg.h"

#include <cmath>


mixtureDistribution::mixtureDistribution(const vector<generalGammaDistribution*>& components, const Vdouble& componentsProb, quadratureType gammaType)
{
	if (components.size() < 1)
		errorMsg::reportError("the number of Gamma components must be positive");
	
	_components.clear();
	for (int i = 0; i < components.size(); ++i)
	{
		generalGammaDistribution* comp = static_cast<generalGammaDistribution*>(components[i]->clone());
		_components.push_back(comp);
	}

	_globalRate = 1.0;
	setComponentsProb(componentsProb);
}


//init the mixture with componentsNum components - the alpha, beta, and probability for each component is assigned "randomly" 
mixtureDistribution::mixtureDistribution(int componentsNum, int categoriesNumInComponent, quadratureType gammaType/*=LAGUERRE*/, MDOUBLE maxAlpha/*=5.0*/, MDOUBLE maxBeta/*=5.0*/)
{
	if (componentsNum < 1)
		errorMsg::reportError("the number of Gamma components must be positive");

	_components.clear();
	Vdouble componentsProb(componentsNum, 0);
	for (int i = 0; i < componentsNum; ++i)
	{
		MDOUBLE alpha = talRandom::giveRandomNumberBetweenZeroAndEntry(maxAlpha);
		MDOUBLE beta = talRandom::giveRandomNumberBetweenZeroAndEntry(maxBeta);
		componentsProb[i] = talRandom::giveRandomNumberBetweenZeroAndEntry(1.0);
		generalGammaDistribution* pComp; 
		switch (gammaType)
		{
		case LAGUERRE:
			pComp = new generalGammaDistributionLaguerre(alpha, beta, categoriesNumInComponent);
			break;
		case QUANTILE:
			pComp = new generalGammaDistribution(alpha, beta, categoriesNumInComponent);
			break;
		default:
			errorMsg::reportError("unknown quadrature type in mixtureDistribution");
		}
		_components.push_back(pComp);
	}

	scaleVec(componentsProb, 1.0/componentsNum);
	setComponentsProb(componentsProb);
	_globalRate = 1.0;
}
//init the mixture with componentsNum components - the alpha, beta, and probability for each component is assigned with given values
mixtureDistribution::mixtureDistribution(int componentsNum, int categoriesNumInComponent,Vdouble AlphaInit ,Vdouble BetaInit, Vdouble componentProbInit ,quadratureType gammaType/*=LAGUERRE*/, MDOUBLE maxAlpha/*=5.0*/, MDOUBLE maxBeta/*=5.0*/)
{
	if (componentsNum < 1)
		errorMsg::reportError("the number of Gamma components must be positive");

	_components.clear();
	Vdouble componentsProb(componentsNum, 0);
	for (int i = 0; i < componentsNum; ++i)
	{
		MDOUBLE alpha = AlphaInit[i];
		MDOUBLE beta = BetaInit[i];
		componentsProb[i] = componentProbInit[i];
		generalGammaDistribution* pComp; 
		switch (gammaType)
		{
		case LAGUERRE:
			pComp = new generalGammaDistributionLaguerre(alpha, beta, categoriesNumInComponent);
			break;
		case QUANTILE:
			pComp = new generalGammaDistribution(alpha, beta, categoriesNumInComponent);
			break;
		default:
			errorMsg::reportError("unknown quadrature type in mixtureDistribution");
		}
		_components.push_back(pComp);
	}

	scaleVec(componentsProb, 1.0/componentsNum);
	setComponentsProb(componentsProb);
	_globalRate = 1.0;
}

mixtureDistribution::mixtureDistribution(const mixtureDistribution& other)
:	_componentsWeight(other._componentsWeight),
	_globalRate(other._globalRate),
	_totalWeight(other._totalWeight)
{
	_components.clear();
	for (int i = 0; i < other.getComponentsNum(); ++i)
	{
		generalGammaDistribution* comp = static_cast<generalGammaDistribution*>(other._components[i]->clone());
		_components.push_back(comp);
	}
}


mixtureDistribution& mixtureDistribution::operator=(const mixtureDistribution &otherDist) 
{
	_globalRate = otherDist._globalRate;
	_componentsWeight = otherDist._componentsWeight;
	_totalWeight = otherDist._totalWeight;
	if (this != &otherDist) // Check for self-assignment
	{
		for (int i = 0; i < getComponentsNum(); ++i)
		{
			if (_components[i] != NULL)
			{
				generalGammaDistribution* pComp = static_cast<generalGammaDistribution*>(otherDist.getComponent(i)->clone());
				delete _components[i];
				_components[i] = pComp;;
			}
		}
	}
	return *this;
}

  
void mixtureDistribution::clear()
{
	for (int i = 0; i < getComponentsNum(); ++i)
	{
		if (_components[i] != NULL)
		{
			delete _components[i];
			_components[i] = NULL;
		}
	}
	_components.clear();
}


mixtureDistribution::~mixtureDistribution()
{
	clear();
}

const int mixtureDistribution::categories() const
{
	int res = 0;
	for (int i = 0; i < getComponentsNum(); ++i)
	{
		res += _components[i]->categories();
	}
	return res;
}

void mixtureDistribution::setComponentsProb(const Vdouble& componentsProb)
{
	if (getComponentsNum() != componentsProb.size())
		errorMsg::reportError("the number of Gamma components is not the same as the number of probabilities");
	_totalWeight = 0.0;
	for (int i = 0; i < componentsProb.size(); ++i)
		_totalWeight += componentsProb[i];
	if (!DEQUAL(_totalWeight, 1.0))
		errorMsg::reportError("the sum of components probabilities must sum to 1.0");
	_componentsWeight = componentsProb;
}


void mixtureDistribution::change_number_of_categoriesPerComp(int in_number_of_categories)
{
	for (int i = 0; i <getComponentsNum(); ++i)
		_components[i]->change_number_of_categories(in_number_of_categories);
}

//change_number_of_components: if the newCompNum is getComponentsNum()-1
//then duplicate one of the components and adjust the probabilities
void mixtureDistribution::change_number_of_components(const int in_number_of_components)
{
	if (getComponentsNum() == in_number_of_components)
		return;
	else if (getComponentsNum() == in_number_of_components - 1)
	{
		//duplicate the first component
		normalizeProbabilities();
		generalGammaDistribution* comp = static_cast<generalGammaDistribution*>(_components[0]->clone());
		_components.push_back(comp);
		//adjust the components probabilities so that the probs of the 
		//two identical components (i.e., 0 and the new Comp) are equal
		_componentsWeight[0] /= 2; 
		_componentsWeight.push_back(_componentsWeight[0]);
		normalizeProbabilities();
	}
	else
		errorMsg::reportError("cannot change the number of components in mixtureDistribution::change_number_of_components()");
}


const MDOUBLE mixtureDistribution::getCumulativeProb(const MDOUBLE x) const
{
	MDOUBLE res = 0.0;
	for (int i = 0; i < getComponentsNum(); ++i)
		res += _components[i]->getCumulativeProb(x) * getComponentProb(i);
	return res;
}

const MDOUBLE mixtureDistribution::rates(const int category) const
{
	if (category > categories() - 1)
		errorMsg::reportError("the required category does not exist!");
	int componentNum, categoryInComponent, totalCat = 0;
	for (int i = 0; i < getComponentsNum(); ++i)
	{
		if (category < _components[i]->categories() + totalCat)
		{
			componentNum = i;
			categoryInComponent = category - totalCat;
			break;
		}
		totalCat += _components[i]->categories();
	}
	return _components[componentNum]->rates(categoryInComponent) * _globalRate;
}

const MDOUBLE mixtureDistribution::ratesProb(const int category) const
{
	if (category > categories() - 1)
		errorMsg::reportError("there required category does not exist!");
	int componentNum, categoryInComponent, totalCat = 0;
	for (int i = 0; i < getComponentsNum(); ++i)
	{
		if (category < _components[i]->categories() + totalCat)
		{
			componentNum = i;
			categoryInComponent = category - totalCat;
			break;
		}
		totalCat += _components[i]->categories();
	}	
	return  getComponentProb(componentNum) * _components[componentNum]->ratesProb(categoryInComponent);
}


void mixtureDistribution::setMixtureParameters(const Vdouble& alphaVec, const Vdouble& betaVec, const Vdouble& componentsProb)
{
	if (alphaVec.size() != getComponentsNum())
		errorMsg::reportError("the size of the alphas vector is not identical to the number of components");
	if (betaVec.size() != getComponentsNum())
		errorMsg::reportError("the size of the batas vector is not identical to the number of components");
	if (componentsProb.size() != getComponentsNum())
		errorMsg::reportError("the size of the components probabilities vector is not identical to the number of components");

	setComponentsProb(componentsProb);
	int categoriesInComponent = _components[0]->categories(); 
	for (int i = 0; i < getComponentsNum(); ++i)
		_components[i]->setGammaParameters(categoriesInComponent, alphaVec[i], betaVec[i]);
}

//the following functions set the components probabilities.
//Note, that the new prob is not inWeight, but is scaled so that the total probabilities are 1.0
void mixtureDistribution::setComponentWeight(MDOUBLE inWeight, const int componentNum, const MDOUBLE minWeight/*=0.01*/)
{
	if((inWeight<0.0) || (inWeight>1.0)){
		errorMsg::reportError("the probability assignment is not [0,1]");
	}
	if (inWeight < minWeight)
		inWeight = minWeight;
	MDOUBLE otherProbs = 1-inWeight;
	Vdouble probs(getComponentsNum(), 0.0);
	MDOUBLE sumOther = 0.0;
	int i;
	for (i = 0; i < getComponentsNum(); ++i)
	{
		if (i != componentNum)
			sumOther += _componentsWeight[i];
	}
	MDOUBLE factor = otherProbs / sumOther;
	for (i = 0; i < getComponentsNum(); ++i)
	{
		probs[i] = _componentsWeight[i] * factor ;
	}
	probs[componentNum] = inWeight;
	setComponentsProb(probs);
 	
	//_totalWeight -= _componentsWeight[componentNum];
 //   _componentsWeight[componentNum] = inWeight;
	//_totalWeight += _componentsWeight[componentNum];
}

//scale the components weights so that they sum to 1.0.
void mixtureDistribution::normalizeProbabilities()
{
	if (_componentsWeight.size() != getComponentsNum())
		errorMsg::reportError("problem in mixtureDistribution::normalizeProbabilities()");
	int i;
	for(i = 0; i < getComponentsNum(); ++i)
	{
		_componentsWeight[i] /= _totalWeight;
	}
	_totalWeight = 1.0;
}

void mixtureDistribution::printParams(ostream& outF)
{		
	MDOUBLE avgRate = 0.0;
	for (int k = 0; k < getComponentsNum(); ++k)
	{
		outF << "comp="<<k<<" Alp/Beta= "<<getAlpha(k)/getBeta(k)<<" alpha= "<<getAlpha(k) << " beta= " <<getBeta(k)<<" Prob= "<<getComponentProb(k)<<endl;  
		avgRate +=  (getAlpha(k) / getBeta(k)) * getComponentProb(k);
	}
	outF<<"# The prior average rate is: " <<avgRate<<endl;
}