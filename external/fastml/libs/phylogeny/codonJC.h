// $Id: codonJC.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___CODON_JC
#define ___CODON_JC

#include "replacementModel.h"
#include <cmath>
using namespace std;

namespace codonDef {
	const MDOUBLE Alp = 61.0;
	const MDOUBLE odAl = 1.0/Alp; // one divided by alphabet
	const MDOUBLE om_odAl = 1.0-odAl; // one minus odAl;
	const MDOUBLE alDiv_omalp = Alp/(Alp-1.0);
	const MDOUBLE m_alDiv_omalp = -alDiv_omalp;
}

class codonJC : public replacementModel {
public:

	virtual replacementModel* clone() const { return new codonJC(*this); }// see note down:
	const int alphabetSize() const {return 61;}

	explicit codonJC(){};
	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const {
		return ((i==j) ?  codonDef::odAl+codonDef::om_odAl*exp(codonDef::m_alDiv_omalp*d): codonDef::odAl-codonDef::odAl*exp(codonDef::m_alDiv_omalp*d));
	}

	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const{
			return ((i==j) ?  -exp(codonDef::m_alDiv_omalp*d): exp(codonDef::m_alDiv_omalp*d)/(codonDef::Alp-1));
	}
	const MDOUBLE freq(const int i) const {return codonDef::odAl;};

	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const{
		return ((i==j) ?  codonDef::alDiv_omalp*exp(codonDef::m_alDiv_omalp*d): codonDef::m_alDiv_omalp*exp(codonDef::m_alDiv_omalp*d));
	}

};

#endif 

// note: according to the new C++ rules, the clone function should be like this:
//	virtual aaJC* clone() const { return new aaJC(*this); }
// however, not all compiler support it yet. look at More Effective C++ page 126.



