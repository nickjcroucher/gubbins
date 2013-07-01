// $Id: tamura92.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___TAMURA92
#define ___TAMURA92

#include "replacementModel.h"
#include <cmath>

class tamura92 : public replacementModel {
public:
	explicit tamura92(const MDOUBLE theta,
					  const MDOUBLE TrTv);

	virtual replacementModel* clone() const { return new tamura92 (*this); }

	const int alphabetSize() const {return 4;}
	inline void changeTrTv(const MDOUBLE TrTv) { _TrTv = TrTv; }
	void changeTheta(const MDOUBLE theta);
	MDOUBLE getTrTv() const {return _TrTv;}
	MDOUBLE getTheta() const {return _theta;}

	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE freq(const int i) const {return _freq[i];};
	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const;

	const MDOUBLE dPij_tdBeta(const int i, const int j, const MDOUBLE t) const;

private:
	Vdouble _freq;
	MDOUBLE _theta;
	MDOUBLE _TrTv;
};

#endif

