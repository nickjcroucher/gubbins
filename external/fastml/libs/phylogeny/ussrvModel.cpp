// 	$Id: ussrvModel.cpp 962 2006-11-07 15:13:34Z privmane $	
#include "ussrvModel.h"

ussrvModel::ussrvModel(const stochasticProcess& baseSp, const stochasticProcessSSRV& ssrvSp, const MDOUBLE& f)
: _f(f),_baseSp(NULL),_ssrvSp(NULL)
{
	_baseSp = new stochasticProcess(baseSp);
	_ssrvSp = new stochasticProcessSSRV(ssrvSp);
	
	// get alpha from sp
	replacementModelSSRV* pMulRM = static_cast<replacementModelSSRV*>(_ssrvSp->getPijAccelerator()->getReplacementModel());
	_alpha = static_cast<gammaDistribution*>(pMulRM->getDistribution())->getAlpha(); 

	// check that alpha is equal the baseSp alpha
	MDOUBLE baseSpAlpha = static_cast<gammaDistribution*>(baseSp.distr())->getAlpha();
	if (_alpha != baseSpAlpha)
		errorMsg::reportError("Error in the constructor of ussrvModel. alpha of the ssrv stochastic process is different from that of the base model");
}

ussrvModel::~ussrvModel()
{
	if (_baseSp) delete _baseSp;
	if (_ssrvSp) delete _ssrvSp;
}

ussrvModel::ussrvModel(const ussrvModel& other)
{
	_f = other._f;
	_baseSp = new stochasticProcess(*other._baseSp);
	_ssrvSp = new stochasticProcessSSRV(*other._ssrvSp);
}

ussrvModel& ussrvModel::operator=(const ussrvModel& other)
{
	if (_baseSp) delete _baseSp;
	if (_ssrvSp) delete _ssrvSp;
	
	_f = other._f;
	_alpha = other._alpha;

	_baseSp = new stochasticProcess(*other._baseSp);
	_ssrvSp = new stochasticProcessSSRV(*other._ssrvSp);
	
	return *this;
}

void ussrvModel::updateAlpha(const MDOUBLE& alpha)
{
	_alpha = alpha;
	if (alpha<0) 
	{
		LOG(4, << "ussrvModel::updateAlpha , alpha is < 0 " << endl);
		return;
	}
	// update alpha of the ssrv model
	replacementModelSSRV* pMulRM = static_cast<replacementModelSSRV*>(_ssrvSp->getPijAccelerator()->getReplacementModel());
	gammaDistribution* gammaDist = static_cast<gammaDistribution*>(pMulRM->getDistribution()); 
	gammaDist->setAlpha(alpha);
	pMulRM->updateQ();

	// update alpha of the base model
	(static_cast<gammaDistribution*>(_baseSp->distr()))->setAlpha(alpha);
}

void ussrvModel::updateNu(const MDOUBLE& nu)
{
	if (nu<0) 
	{
		LOG(4,<<"ussrvModel::updateNu , nu is < 0 " <<endl);
		return;
	}
	static_cast<replacementModelSSRV*>(_ssrvSp->getPijAccelerator()->getReplacementModel())->setRateOfRate(nu);
}

MDOUBLE ussrvModel::getNu() const 
{
	return (static_cast<replacementModelSSRV*>(_ssrvSp->getPijAccelerator()->getReplacementModel())->getRateOfRate());
}

void ussrvModel::updateF(const MDOUBLE& f) 
{
	if ((f<0) || (f>1))
	{
		LOG(4,<<"ussrvModel::updateF , f must be between 0 to 1. f is: "<< f << endl);
		return;
	}
	_f=f;
}

// In order for the branch lengths and the nu parameter to be meaningfull, one must normalize the
// matrices of both the replacement models (the base model and the ssrv model)
// so that f*Sigma[i](PiQij) + (1-f)*Sigma[i](P`iQ`ij) = 1 (for i!=j)
// where Q and P belong to the ssrv model, P` and Q` belong to the base model. (Q` doesn't include the rates)
// The normalization doesn't affect the likelihood.
// see below for more explanations.
// Theoretically, we should therefore calculate this weighted sumPijQij (Denote by x), and then:
// 1) devide nu by x.
// 2) devide all the rates (of the base model and of the ssrv model) by x. 
// (this could be done using the _globalRate member of the gammaDistribution class)
// 3) multiply every branch length by x.
// Instead, we just report x, so that the user can do all this whenever he wishes to.

MDOUBLE ussrvModel::calcNormalizeFactor()
{
	// calculate sumPijQij
	MDOUBLE sumPijQij = 0.0;
	int i;
	// of the base model
	int baseAlphabetSize = _baseSp->alphabetSize();
	for (i=0; i < baseAlphabetSize; ++i)
		sumPijQij-= _baseSp->freq(i) * _baseSp->dPij_dt(i,i,0);
	sumPijQij*=(1-_f);
	
	// of the ssrv model
	sumPijQij+=_f*static_cast<replacementModelSSRV*>(_ssrvSp->getPijAccelerator()->getReplacementModel())->sumPijQij();	

	return sumPijQij;
}

// This is not done when using normal sp (instead of ussrvModel), since:
// average(rates)=1 --> 
// (for 2 categories, f=0.5, 1-f =0.5) 0.5*r1*Sigma[i](PiQij) + 0.5*r2*Sigma[i](PiQij) = 1 --> 
// (since (r1+r2)*0.5 = 1) Sigma[i](PiQij) = 1 . This is always true, and taken care of in the readMatrix
// method.

