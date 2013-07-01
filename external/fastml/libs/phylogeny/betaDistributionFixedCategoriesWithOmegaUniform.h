#ifndef ___BETA_DISTR_FIXED_CATEGORIES_OMEGA_UNIFORM
#define ___BETA_DISTR_FIXED_CATEGORIES_OMEGA_UNIFORM
/************************************************************
This class differ from the regular betaOmegaDistribution in that 
the rateCategories are fixed according to the user's decision. 
Thus, only the probability of each category changes for each specific alpha value but 
the rate categories themselves are constant.
************************************************************/
#include "definitions.h"
#include "betaDistributionFixedCategories.h"
#include "uniformDistribution.h"
#include "errorMsg.h"


class betaDistributionFixedCategoriesOmegaUniform : public distribution {
public:
	
	explicit betaDistributionFixedCategoriesOmegaUniform(const betaDistributionFixedCategoriesOmegaUniform& other);
	explicit betaDistributionFixedCategoriesOmegaUniform(int betaDistrCatNum,MDOUBLE alpha,MDOUBLE beta,
														int omegaCatNum =10,MDOUBLE omegaLowerBound = 1,MDOUBLE omegaUpperBound = 11);
	explicit betaDistributionFixedCategoriesOmegaUniform() {};
	virtual ~betaDistributionFixedCategoriesOmegaUniform() {};
	virtual distribution* clone() const { return new betaDistributionFixedCategoriesOmegaUniform(*this); }
	virtual void change_number_of_categories(int in_number_of_categories);
	virtual void setBetaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta);
	
	virtual const int categories() const {return _betaDistr.categories()+ _omegaDistr.categories();}
	virtual const int betaCategories()const {return _betaDistr.categories();};
	virtual const MDOUBLE rates(const int i) const;
	virtual const MDOUBLE ratesProb(const int i_rate) const;
	virtual void setGlobalRate(const MDOUBLE x) {_betaDistr.setGlobalRate(x);}
	virtual MDOUBLE getGlobalRate()const {return  _betaDistr.getGlobalRate();}
	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const;
	virtual void setAlpha(MDOUBLE newAlpha){ _betaDistr.setAlpha(newAlpha);}
	virtual MDOUBLE getAlpha() const {return _betaDistr.getAlpha();};
	virtual void setBeta(MDOUBLE newBeta){_betaDistr.setBeta(newBeta);}
	virtual MDOUBLE getBeta() const {return _betaDistr.getBeta();};
	virtual MDOUBLE getBorder(const int i) const {return _betaDistr.getBorder(i);}	//return the ith border. Note:  _bonderi[0] = 0, _bondery[categories()] = infinite
	//virtual MDOUBLE getOmegai() const ;
	//virtual MDOUBLE getBetaProbi() const ;
	//virtual void setOmegai(MDOUBLE omega);
	//virtual void setBetaProbi(MDOUBLE betaProb);


private:	
	betaDistributionFixedCategories _betaDistr; //10 fixed cat 0.05, 0.15, 0.25 ...,0.95 
	uniformDistribution _omegaDistr; // w ~ U(1,11) with 10 cat
};



#endif

