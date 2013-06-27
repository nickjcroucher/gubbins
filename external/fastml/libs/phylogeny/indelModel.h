// 	$Id: indelModel.h 962 2006-11-07 15:13:34Z privmane $	
#ifndef ___INDEL_MODEL
#define ___INDEL_MODEL

#include "replacementModel.h"
#include <cmath>
using namespace std;

class indelModel : public replacementModel
{
public:
	explicit indelModel(const MDOUBLE freq_x, const MDOUBLE freq_g)
	{
		_alpha = 1/(2*freq_x*freq_g);
		_freq.push_back(freq_x);
		_freq.push_back(freq_g);
	}
	
	virtual const MDOUBLE Pij_t(const int i, const int j, const MDOUBLE t) const
	{
		if (i==j)
			return exp(-t*_alpha);
		return (1-exp(-t*_alpha));
	}
	
	virtual const MDOUBLE freq(const int i) const { return _freq[i];}
	
	virtual const MDOUBLE dPij_dt(const int i, const int j, const MDOUBLE t) const
	{
		// [e^(-t/2PxPg)] / 2PxPg
		return (exp(-t*_alpha)*_alpha);
	}
	virtual const MDOUBLE d2Pij_dt2(const int i, const int j, const MDOUBLE t) const
	{
		// [-e^(-t/2PxPg)] / [(2PxPg)^2]
		return ( -exp(-t*_alpha) * _alpha * _alpha);
	}

	virtual replacementModel* clone() const { return new indelModel(*this);}

	virtual	const int alphabetSize() const {return 2;};


	void setFreqX(const MDOUBLE freq_x); 
	void setFreqG(const MDOUBLE freq_g); 


private:
	Vdouble _freq; // [0] X [1] -  
	// save _alpha to make things faster. _alpha depends on _freq
	MDOUBLE _alpha;
};


#endif 


 
	
	

