
#ifndef ___DIST_PLUS_CATEGORY
#define ___DIST_PLUS_CATEGORY

#include "definitions.h"
#include "distribution.h"
#include "logFile.h"
#include "errorMsg.h"

class distributionPlusCategory : public distribution {

public:
	explicit distributionPlusCategory(const distribution* pBaseDist, MDOUBLE baseDistProb,MDOUBLE categoryVal,MDOUBLE globalRate=1);
	explicit distributionPlusCategory();
	explicit distributionPlusCategory(const distributionPlusCategory& other);
	virtual ~distributionPlusCategory();
	virtual distributionPlusCategory& operator=(const distributionPlusCategory &other);
	virtual distribution* clone() const { return new distributionPlusCategory(*this); }

	distribution* getBaseDistribution() {return _pBaseDist;}
	virtual const int categories() const;
	virtual const MDOUBLE rates(const int category) const;
	virtual const MDOUBLE ratesProb(const int category) const;

	virtual void setGlobalRate(const MDOUBLE x) {_globalRate=x;}
	virtual MDOUBLE getGlobalRate()const {return  _globalRate;}
	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const;
	virtual void change_number_of_categories(int in_number_of_categories);

	virtual MDOUBLE getCategoryVal() const {return _categoryVal;}
	virtual MDOUBLE getBaseDistProb() const {return _baseDistProb;}
	virtual void setCategoryVal(MDOUBLE categoryVal) { _categoryVal = categoryVal;}
	virtual void setBaseDistProb(MDOUBLE baseDistProb); 

protected:	
	MDOUBLE _globalRate;
	distribution* _pBaseDist;
	MDOUBLE _categoryVal;
	MDOUBLE _baseDistProb;

};

#endif // ___DIST_PLUS_CATEGORY
