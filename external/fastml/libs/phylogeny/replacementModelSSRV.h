// 	$Id: replacementModelSSRV.h 1914 2007-04-04 08:40:35Z osnatz $	
#ifndef ___REPLACEMENT_MODEL_SSRV
#define ___REPLACEMENT_MODEL_SSRV

#include <cmath>
#include "replacementModel.h"
#include "distribution.h"
#include "fromQtoPt.h"
#include "errorMsg.h"
#include "definitions.h"

class replacementModelSSRV : public replacementModel 
{
public:
	explicit replacementModelSSRV(const distribution* dist, const replacementModel* baseRM, MDOUBLE rateOfRate = 1);
	explicit replacementModelSSRV(const replacementModelSSRV& other);
	~replacementModelSSRV();
	replacementModelSSRV& operator=(const replacementModelSSRV &other); 
	const int alphabetSize() const;
	virtual replacementModel* clone() const  {return new replacementModelSSRV(*this);}
	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const {
		return _q2pt.Pij_t(i,j,d);
	}
	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const{
		return _q2pt.dPij_dt(i,j,d);
	}
	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const{
		return _q2pt.d2Pij_dt2(i,j,d);
	}
	
	const MDOUBLE freq(const int i) const {return _freq[i];}
	
	distribution* getDistribution() const { return _dist;} // @@@@ this const is a lie !!!
	void setDistribution(const distribution* dist);  // it's important to call updateQ after changing the distribution parameters
	
	replacementModel* getBaseRM() const { return _baseRM;} // @@@@ this const is a lie (for the same reason as getDistribution()

	MDOUBLE getRateOfRate() const { return _rateOfRate;}
	void setRateOfRate(MDOUBLE rateOfRate) { _rateOfRate=rateOfRate; updateQ();}

	VVdouble getQ() const { return _Q;}
	Vdouble getFreqs() const {return _freq;}

	MDOUBLE sumPijQij() const;

	void updateQ();
	void updateFreq();
	q2pt getQ2pt() const {return _q2pt;} // used for debug only
	//void norm(MDOUBLE scale);

private:
	distribution* _dist;
	replacementModel* _baseRM;
	MDOUBLE _rateOfRate;
	q2pt _q2pt;
	Vdouble _freq;
	VVdouble _Q;
			
};

#endif 

/* @@@@ When we want to optimize alpha, we usually get the distibution from the stochastic process and then
convert it using static_cast, for example to gammaDistribution and use its method setAlpha.
For this reason, the method distr() in replacmentModel and the method getDistribution here are both const, although 
they actually allow changing the distribution.
A good solution for this is to add a setDistribution in the stochasticProcess.
This will check if the distributions are of the same type and if so, will just update the alpha.
*/

// @@@@ Idea - maybe there is no need of replacementModelSSRV. This can be stochasticProcessSSRV - not good. the SP also has an accelerator.


