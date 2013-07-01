// $Id: alphaTrivialAccelerator.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___ALPHA_TRIVIAL_ACCELERATOR
#define ___ALPHA_TRIVIAL_ACCELERATOR

#include "pijAccelerator.h"
#include "readDatMatrix.h"
class alphaTrivialAccelerator : public pijAccelerator {
public:

	explicit alphaTrivialAccelerator(pupAll* pb, const MDOUBLE alpha) :
		_pb(static_cast<pupAll *> (pb->clone())),
		_alpha(alpha) 
		{};

	alphaTrivialAccelerator(const alphaTrivialAccelerator& other):
		_pb(NULL),
		_alpha(other._alpha) {
		if (other._pb != NULL) 
			_pb = static_cast<pupAll *>(other._pb->clone());
	}

	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const {return _pb->Pij_tAlpha(i,j,d,_alpha);}

	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const{return _pb->Pij_tAlpha_dt(i,j,d,_alpha);};

	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const{return _pb->Pij_tAlpha_dt2(i,j,d,_alpha);};

	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d, const MDOUBLE alpha) const {return _pb->Pij_tAlpha(i,j,d,alpha);}

	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d, const MDOUBLE alpha) const{return _pb->Pij_tAlpha_dt(i,j,d,alpha);};

	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d, const MDOUBLE alpha) const{return _pb->Pij_tAlpha_dt2(i,j,d,alpha);};

	const MDOUBLE freq(const int i) const{return _pb->freq(i);}

	virtual pijAccelerator* clone() const { return new alphaTrivialAccelerator(*this);}

	virtual ~alphaTrivialAccelerator() {delete _pb;}

	virtual const int alphabetSize() const {return _pb->alphabetSize();}

	virtual replacementModel* getReplacementModel() const {
		return (static_cast<replacementModel * const>(_pb));
	}

	const MDOUBLE alpha(void) const {return _alpha;}
	void setAlpha(const MDOUBLE alpha) {_alpha=alpha;}

private:
	pupAll* _pb;
	MDOUBLE _alpha;
};

#endif

