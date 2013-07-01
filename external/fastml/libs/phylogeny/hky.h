// $Id: hky.h 4291 2008-06-23 10:23:10Z itaymay $

#ifndef ___HKY
#define ___HKY

#include "replacementModel.h"
#include <cmath>

class hky : public replacementModel {
public:
	explicit hky(const MDOUBLE inProb_a,
					const MDOUBLE inProb_c,
					const MDOUBLE inProb_g,
					const MDOUBLE inProb_t,
					const MDOUBLE TrTv);

	explicit hky(vector<MDOUBLE> inProbs, const MDOUBLE TrTv);

	virtual replacementModel* clone() const { return new hky(*this); }
//	virtual nucJC* clone() const { return new nucJC(*this); } // see note down:

	const int alphabetSize() const {return 4;}


	void changeTrTv(const MDOUBLE In_TrTv);
	MDOUBLE getTrTv() const;
	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE freq(const int i) const {return _freq[i];};
	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const;

	const MDOUBLE dPij_tdBeta(const int i, const int j, const MDOUBLE t) const;

private:
	void initParams(MDOUBLE TrTv); // init _a, _b, _c, and _y by using _freq and TrTv

private:
	Vdouble _freq;
	MDOUBLE _a; //
	MDOUBLE _b; //

	MDOUBLE _c,_y; // relationship between probA, probC, prob G, prob T.
};

#endif

