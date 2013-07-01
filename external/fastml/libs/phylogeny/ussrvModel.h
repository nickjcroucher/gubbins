// 	$Id: ussrvModel.h 962 2006-11-07 15:13:34Z privmane $	
#ifndef _USSRV_MODEL
#define _USSRV_MODEL

#include "stochasticProcessSSRV.h"
#include "stochasticProcess.h"
#include "errorMsg.h"
#include "gammaDistribution.h"
#include "replacementModelSSRV.h"
#include "logFile.h"
class ussrvModel
{
public:
	explicit ussrvModel(){errorMsg::reportError("This constractor shold never be used");}
	explicit ussrvModel(const stochasticProcess& baseSp, const stochasticProcessSSRV& ssrvSp, const MDOUBLE& f);
	virtual ~ussrvModel();
	explicit ussrvModel(const ussrvModel& other);
	ussrvModel& operator=(const ussrvModel& other);
	// const int alphabetSize() const ;
	MDOUBLE getF() const {return _f;}
	MDOUBLE getAlpha() const {return _alpha;}
	MDOUBLE getNu() const ;
	const stochasticProcessSSRV& getSSRVmodel() const {return *_ssrvSp;}
	const stochasticProcess& getBaseModel() const {return *_baseSp;}
	int noOfCategor() const {return _baseSp->categories();}
	MDOUBLE getCategorProb(int i)  const {return _baseSp->distr()->ratesProb(i);}

	void updateF(const MDOUBLE& f);
	void updateAlpha(const MDOUBLE& alpha);
	void updateNu(const MDOUBLE& nu); 
	
	MDOUBLE calcNormalizeFactor(); // return the factor according to which the model should be normalized.

private:
	MDOUBLE _f; //probability of SSRV model. The probability of the base model, i.e. no SSRV, is 1-_f .
	MDOUBLE _alpha; // should be always equal to the _baseSp alpha and the _ssrvSp alpha.
	stochasticProcess* _baseSp; // for the base model
	stochasticProcessSSRV* _ssrvSp; // for the SSRV model
};

#endif // _USSRV_MODEL
