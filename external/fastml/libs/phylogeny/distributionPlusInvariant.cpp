#include "definitions.h"
#include "distributionPlusInvariant.h"
#include "errorMsg.h"
#include "logFile.h"

//#define RATE_INVARIANT 1e-10


distributionPlusInvariant::distributionPlusInvariant(
	distribution* pDist, const MDOUBLE pInv, const MDOUBLE globalRate, MDOUBLE rateInvariantVal) 
{
	_globalRate=globalRate;
	_Pinv = pInv;
	_rateInvariantVal = rateInvariantVal;
	_pBaseDist = NULL;
	if (pDist!= NULL)
		_pBaseDist = pDist->clone();
}

distributionPlusInvariant::distributionPlusInvariant() 
{
	_globalRate=1.0;
	_Pinv = 0;
	_rateInvariantVal = 0;
	_pBaseDist = NULL;
}


distributionPlusInvariant& distributionPlusInvariant::operator=(const distributionPlusInvariant& other)
{
	_globalRate = other._globalRate;
	_Pinv = other._Pinv;
	_rateInvariantVal = other._rateInvariantVal;
	_pBaseDist = NULL;	
	if (other._pBaseDist != NULL) 
		_pBaseDist = other._pBaseDist->clone();
	return *this;
}

distributionPlusInvariant::~distributionPlusInvariant()
{
	if (_pBaseDist != NULL)
		delete _pBaseDist;
}


//gets cumulative probability till a certain point
const MDOUBLE distributionPlusInvariant::getCumulativeProb(const MDOUBLE x) const
{
	if (x < 0)
		errorMsg::reportError("x < 0 in distributionPlusInvariant::getCumulativeProb()");
	return (_Pinv + (1 -_Pinv) * _pBaseDist->getCumulativeProb(x));
}


const MDOUBLE distributionPlusInvariant::ratesProb(const int category) const
{
	if (category == categories()-1)
		return _Pinv;
	else
		return (1 - _Pinv) * _pBaseDist->ratesProb(category);
}

const MDOUBLE distributionPlusInvariant::rates(const int category) const
{
	if (category == categories()-1)
		return _rateInvariantVal; //RATE_INVARIANT
	else
		return _pBaseDist->rates(category);
}

const int distributionPlusInvariant::categories() const
{
	return 1 + _pBaseDist->categories(); 
}


