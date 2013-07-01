// $Id: trivialAccelerator.h 1925 2007-04-04 16:40:22Z privmane $

#ifndef ___TRIVIAL_ACCELERATOR
#define ___TRIVIAL_ACCELERATOR

#include "pijAccelerator.h"
#include "replacementModel.h"

class trivialAccelerator : public pijAccelerator {
public:

	explicit trivialAccelerator(const replacementModel* pb): _pb(pb->clone()) {};
	trivialAccelerator(const trivialAccelerator& other):_pb(NULL){if (other._pb != NULL) _pb = other._pb->clone();}
	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const {return _pb->Pij_t(i,j,d);}
	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const{return _pb->dPij_dt(i,j,d);};
	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const{return _pb->d2Pij_dt2(i,j,d);};
	const MDOUBLE freq(const int i) const{return _pb->freq(i);}
	virtual pijAccelerator* clone() const { return new trivialAccelerator(*this);}
	virtual ~trivialAccelerator() {delete _pb;}
	virtual const int alphabetSize() const {return _pb->alphabetSize();}
	virtual replacementModel* getReplacementModel() const {return (_pb);}

private:
	replacementModel* _pb;
};

#endif

// There is no distribution in the trivial accelerator. Actually, it's just an interface 
// to the replacement Model and it doesn't accelerate anything.
// Every method retruns exactly the replacementModel corresponding method result.

