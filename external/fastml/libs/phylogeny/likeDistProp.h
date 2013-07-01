// $Id: likeDistProp.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___LIKE_DIST_PROP
#define ___LIKE_DIST_PROP

#include "definitions.h"
#include "countTableComponent.h"
#include "stochasticProcess.h"
#include <cmath>

class likeDistProp {
private:
	const int _alphabetSize;
	const vector<stochasticProcess>& _s1;
	const MDOUBLE _toll;
public:
	const MDOUBLE giveDistance(	const vector<countTableComponentGam>& ctc,
								MDOUBLE& resL) const;
	explicit likeDistProp(const int alphabetSize,
							const vector<stochasticProcess>& s1,
							const MDOUBLE toll) : _alphabetSize(alphabetSize), _s1(s1) ,_toll(toll){
	}
};



class C_evallikeDistProp_d{ // derivative.
public:
  C_evallikeDistProp_d(const vector<countTableComponentGam>& ctc,
				 const vector<stochasticProcess>& inS1)    : _ctc(ctc), _sp(inS1) {};
private:
	const  vector<countTableComponentGam>& _ctc;
	const vector<stochasticProcess>& _sp;
public:
	MDOUBLE operator() (MDOUBLE dist) {
		MDOUBLE	sumDL=0.0;
		const MDOUBLE epsilonPIJ = 1e-10;
		for (int gene=0; gene < _ctc.size(); ++ gene) {
			for (int alph1=0; alph1 <  _ctc[gene].alphabetSize(); ++alph1){
				for (int alph2=0; alph2 <  _ctc[gene].alphabetSize(); ++alph2){
					for (int rateCategor = 0; rateCategor<_sp[gene].categories(); ++rateCategor) {
						MDOUBLE rate = _sp[gene].rates(rateCategor);
						MDOUBLE pij= _sp[gene].Pij_t(alph1,alph2,dist*rate);
						MDOUBLE dpij = _sp[gene].dPij_dt(alph1,alph2,dist*rate);
						if (pij<epsilonPIJ) {
							pij = epsilonPIJ;
							dpij = epsilonPIJ;
						}
						sumDL+= _ctc[gene].getCounts(alph1,alph2,rateCategor)*dpij*_sp[gene].ratesProb(rateCategor)
										*rate/pij;
					}
				}
			}
		}
		return -sumDL;
	}
};



class C_evallikeDistProp{
private:
	const vector<countTableComponentGam>& _ctc;
	const vector<stochasticProcess>& _sp;
public:
	C_evallikeDistProp(const vector<countTableComponentGam>& ctc,
					const vector<stochasticProcess>& inS1):_ctc(ctc), _sp(inS1) {};

	MDOUBLE operator() (MDOUBLE dist) {
		const MDOUBLE epsilonPIJ = 1e-10;
		MDOUBLE sumL=0.0;
		for (int gene=0; gene < _ctc.size(); ++ gene) {
			for (int alph1=0; alph1 < _ctc[gene].alphabetSize(); ++alph1){
				for (int alph2=0; alph2 <  _ctc[gene].alphabetSize(); ++alph2){
					for (int rateCategor = 0; rateCategor<_sp[gene].categories(); ++rateCategor) {
						MDOUBLE rate = _sp[gene].rates(rateCategor);
						MDOUBLE pij= _sp[gene].Pij_t(alph1,alph2,dist*rate);
						if (pij<0) {
							pij = epsilonPIJ;
						}
						sumL += _ctc[gene].getCounts(alph1,alph2,rateCategor)*(log(pij)-log(_sp[gene].freq(alph2)))*_sp[gene].ratesProb(rateCategor);
					}
				}
			}
		}
		return -sumL;
	}
};

#endif

