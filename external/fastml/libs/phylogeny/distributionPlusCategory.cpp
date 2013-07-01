#include "distributionPlusCategory.h"

distributionPlusCategory::distributionPlusCategory(const distribution* pBaseDist, MDOUBLE baseDistProb,MDOUBLE categoryVal,MDOUBLE globalRate)
:
_globalRate(globalRate),
_categoryVal(categoryVal),
_baseDistProb(baseDistProb)
{
	if (pBaseDist!= NULL)
		_pBaseDist = pBaseDist->clone();
}

distributionPlusCategory::distributionPlusCategory()
:
_globalRate(1.0),
_pBaseDist(NULL),
_categoryVal(1.0),
_baseDistProb(0.0)
{
}

distributionPlusCategory::distributionPlusCategory(const distributionPlusCategory& other)
{
	(*this) = other;
}

distributionPlusCategory& distributionPlusCategory::operator=(const distributionPlusCategory &other)
{
	if (this != &other) 
	{
		_globalRate = other._globalRate;
		if (other._pBaseDist) {
			_pBaseDist = other._pBaseDist->clone();
		}
		else {
			_pBaseDist = NULL;
		}
		_categoryVal = other._categoryVal;
		_baseDistProb = other._baseDistProb;

	}
   return *this;
}

distributionPlusCategory::~distributionPlusCategory()
{
	if (_pBaseDist)
		delete _pBaseDist;
}

const int distributionPlusCategory::categories() const
{
	 return _pBaseDist->categories()+1;
}


const MDOUBLE distributionPlusCategory::rates(const int category) const
{
	if (category < _pBaseDist->categories())
		return _pBaseDist->rates(category);
	else 
		return _categoryVal;
}


const MDOUBLE distributionPlusCategory::ratesProb(const int category) const
{
	if (category < _pBaseDist->categories())
		return _pBaseDist->ratesProb(category) * _baseDistProb;
	else 
		return (1-_baseDistProb); //category prob
}


//gets cumulative probability till a certain point
const MDOUBLE distributionPlusCategory::getCumulativeProb(const MDOUBLE x) const
{
	MDOUBLE res(0.0);
	if (x < 0)
		errorMsg::reportError("x < 0 in distributionPlusCategory::getCumulativeProb()");
	if (x > _categoryVal - EPSILON)
		res += 1-_baseDistProb;
	res += _baseDistProb * _pBaseDist->getCumulativeProb(x);
	return res;
}


void distributionPlusCategory::change_number_of_categories(int in_number_of_categories)
{
	_pBaseDist->change_number_of_categories(in_number_of_categories);
}


void distributionPlusCategory::setBaseDistProb(MDOUBLE baseDistProb)
{
	if ((baseDistProb < 0.0) || (baseDistProb>1.0) )
		errorMsg::reportError("illegal baseDistProb in distributionPlusCategory::setBaseDistProb");
	
	_baseDistProb = baseDistProb;
}