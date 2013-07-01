#ifndef ___MIXTURE_DIST
#define ___MIXTURE_DIST
/************************************************************
The mixture distribution is combined of several gamma distributions (components).
Each one of the gamma component has its own probability of occurance = Hi, 
such that the sum of Hi equals 1.0.
The categories probabilities are the probability of each component multiply by the category probabilty in the component.
In case the Laguerre option is on: 
the actuall number of cateories (per component) can be lower than the requested number of categories.
************************************************************/
#include "definitions.h"
#include "generalGammaDistribution.h"

class mixtureDistribution : public distribution {
public:
	explicit mixtureDistribution(const vector<generalGammaDistribution*>& components, const Vdouble& componentsProb, quadratureType gammaType);
	explicit mixtureDistribution(int componentsNum, int categoriesNumInComponent, quadratureType gammaType = LAGUERRE, MDOUBLE maxAlpha = 15.0, MDOUBLE maxBeta = 15.0);
	explicit mixtureDistribution(int componentsNum, int categoriesNumInComponent,Vdouble AlphaInit ,Vdouble BetaInit, Vdouble componentProbInit ,quadratureType gammaType = QUANTILE, MDOUBLE maxAlpha = 15.0, MDOUBLE maxBeta = 15.0);

	mixtureDistribution(const mixtureDistribution& other);	

	mixtureDistribution& operator=(const mixtureDistribution &otherDist); 
	virtual distribution* clone() const { return new mixtureDistribution(*this); }
	virtual ~mixtureDistribution();

	//get+set the parameters of the mixture
	void setMixtureParameters(const Vdouble& alphaVec, const Vdouble& betaVec, const Vdouble& componentsProb);
	const generalGammaDistribution* getComponent(int componentNum) const {return _components[componentNum];}
	const int getComponentsNum() const  {return _components.size();}
	const int categories() const; 
	//change_number_of_categoriesPerComp: change the number of categorites for each component. The total number of categories will be (in_number_of_categories*componentNum)
	void change_number_of_categoriesPerComp(int in_number_of_categories);
	void change_number_of_components(const int in_number_of_components);
	const int categoriesForOneComponent() const {return _components[0]->categories();}
	MDOUBLE getAlpha(int componentNum) const {return _components[componentNum]->getAlpha();}
	void setAlpha(MDOUBLE newAlpha, int componentNum) {_components[componentNum]->setAlpha(newAlpha);}
	MDOUBLE getBeta(int componentNum) const {return _components[componentNum]->getBeta();}
	void setBeta(MDOUBLE newBeta, int componentNum) {_components[componentNum]->setBeta(newBeta);}
	void setGammaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta, int componentNum) {_components[componentNum]->setGammaParameters(numOfCategories ,alpha, beta);}
	const MDOUBLE getComponentProb(int componentNum) const {return _componentsWeight[componentNum] / _totalWeight;} 
	void setComponentsProb(const Vdouble& componentsProb);
 	void setGlobalRate(const MDOUBLE r) {_globalRate = r;}
 	MDOUBLE getGlobalRate() const {return _globalRate;}

	//the following function set the components weights.
	//Note that the new component prob is not inWeight, but is scaled so that the total probabilities are 1.0
	void setComponentWeight(MDOUBLE inWeight, const int componentNum, const MDOUBLE minWeight =0.01);
	const MDOUBLE getComponentWeight(int componentNum) const {return _componentsWeight[componentNum];} 
	//scale the components weights so that they sum to 1.0.
	void normalizeProbabilities();

	//get distribution statistics
	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const;
	virtual const MDOUBLE rates(const int category) const;
	virtual const MDOUBLE ratesProb(const int i) const;

	void printParams(ostream& outF );

private:
	void clear();
private:	
	vector<generalGammaDistribution*> _components;
	Vdouble _componentsWeight;
	MDOUBLE _globalRate;
	MDOUBLE _totalWeight; //holds the sum of the components probabilities. This is saved so that we don't need to sum all weight each time getProb() is called
};
#endif
