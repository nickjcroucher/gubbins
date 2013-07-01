// $Id: nucJC.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___NUC_JC
#define ___NUC_JC

#include <cmath>
#include "replacementModel.h"

namespace nucDef {
	const MDOUBLE Alp = 4.0;
	const MDOUBLE odAl = 1.0/Alp; // one divided by alphabet
	const MDOUBLE om_odAl = 1.0-odAl; // one minus odAl;
	const MDOUBLE alDiv_omalp = Alp/(Alp-1.0);
	const MDOUBLE m_alDiv_omalp = -alDiv_omalp;
}

class nucJC : public replacementModel {
public:
	const int alphabetSize() const {return 4;}

	virtual replacementModel* clone() const { return new nucJC(*this); }

	explicit nucJC(){};
	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const {
//		return ((i==j) ?  0.25+0.75*exp(-4.0/3.0*d): 0.25-0.25*exp(-4.0/3.0*d));
		return ((i==j) ?  nucDef::odAl+nucDef::om_odAl*exp(nucDef::m_alDiv_omalp*d): nucDef::odAl-nucDef::odAl*exp(nucDef::m_alDiv_omalp*d));
	}

	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const{
//		return ((i==j) ?  -exp(-4.0/3.0*d): exp(-4.0/3.0*d)/3.0);
		return ((i==j) ?  -exp(nucDef::m_alDiv_omalp*d): exp(nucDef::m_alDiv_omalp*d)/(nucDef::Alp-1));
	}
	const MDOUBLE freq(const int i) const {return 0.25;};

	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const{
	//	return ((i==j) ?  4.0/3.0*exp(-4.0/3.0*d): -4.0/3.0*exp(-4.0/3.0*d));
		return ((i==j) ?  nucDef::alDiv_omalp*exp(nucDef::m_alDiv_omalp*d): nucDef::m_alDiv_omalp*exp(nucDef::m_alDiv_omalp*d));
	}
	
    const MDOUBLE Q(const int i, const int j) const {
		return ((i == j) ? ( - 1.0) : (1.0 / 3.0));
	}


};

#endif 

// note: according to the new C++ rules, the clone function should be like this:
//	virtual nucJC* clone() const { return new nucJC(*this); }
// however, not all compiler support it yet. look at More Effective C++ page 126.


