// $Id: replacementModel.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___REPLACEMENT_MODEL
#define ___REPLACEMENT_MODEL

#include "definitions.h"

class replacementModel{
public:
	virtual const MDOUBLE Pij_t(const int i, const int j, const MDOUBLE t) const = 0;
	virtual const MDOUBLE freq(const int i) const = 0;
	virtual const MDOUBLE dPij_dt(const int i, const int j, const MDOUBLE t) const =0;
	virtual const MDOUBLE d2Pij_dt2(const int i, const int j, const MDOUBLE t) const =0;
	virtual replacementModel* clone() const = 0;
	virtual ~replacementModel()=0;
	virtual	const int alphabetSize() const =0;

    //virtual const MDOUBLE Q(const int i, const int j, const MDOUBLE r = 1.0) const = 0;
	//note that we ask that sigma over i sigma over j!=i of p(i)Qij = 1;
	//this is beacuse we ask the [sigma over i sigma over j!=i p(i)*pij(d)]/d approaches
	//1 as d -> 0. (and l'hopital from here).
};


#endif 

