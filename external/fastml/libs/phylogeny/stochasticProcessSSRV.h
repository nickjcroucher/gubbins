// 	$Id: stochasticProcessSSRV.h 1923 2007-04-04 16:38:14Z privmane $	


#ifndef ___STOCHASTIC_PROCESS_SSRV
#define ___STOCHASTIC_PROCESS_SSRV

#include "stochasticProcess.h"
#include "replacementModelSSRV.h"

// This is a Stochastic process that its distribution is located inside its accelerator.
// _dist should be NULL all the time.
// The number of categories is always 1.
// _pijAccelerator must contain a replacementModelSSRV* as a member.
// The distribution is located inside the replacement model which is a member of _pijAccelerator.

class stochasticProcessSSRV : public stochasticProcess{
public:
	explicit stochasticProcessSSRV(const pijAccelerator *pijAccelerator) : 
		stochasticProcess() { _pijAccelerator = pijAccelerator->clone();}
	explicit stochasticProcessSSRV() : stochasticProcess() {}
	stochasticProcessSSRV(const stochasticProcessSSRV& other) : stochasticProcess(other) {}
	stochasticProcessSSRV& operator=(const stochasticProcessSSRV &other) {stochasticProcess::operator=(other); return (*this);}
	virtual stochasticProcess* clone() const {return new stochasticProcessSSRV(*this);}
	
	virtual ~stochasticProcessSSRV() {}
	
	virtual const int categories() const { return 1; }
	virtual const MDOUBLE rates(const int i) const {return 1.0;}  
	virtual const MDOUBLE ratesProb(const int i) const {return 1.0;} 

	virtual const MDOUBLE Pij_t(const int i, const int j, const MDOUBLE t) const {
		// as opposed to normal stochastic-process. even when t=0 and i!=j the result might be > 0
		return _pijAccelerator->Pij_t(i,j,t); 
	}
	
	virtual distribution* distr() const; // @@@@ this const is a lie !!!
	virtual void setDistribution(const distribution* in_distr);
	
 	virtual void setGlobalRate(const MDOUBLE x) {distr()->setGlobalRate(x);} // @@@@ should this also call updateQ of the RM ??? Doesn't really metter when using gamma distribution
 	virtual MDOUBLE getGlobalRate() const {return distr()->getGlobalRate();}

	void setRateOfRate(MDOUBLE rateOfRate) { 
		static_cast<replacementModelSSRV*>(_pijAccelerator->getReplacementModel())
			->setRateOfRate(rateOfRate);
	}
};

#endif
