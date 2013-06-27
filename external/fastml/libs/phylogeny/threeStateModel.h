#ifndef ___3STATE_MODEL
#define ___3STATE_MODEL

#include "definitions.h"
#include "replacementModel.h"
#include "fromQtoPt.h"
#include "errorMsg.h"
#include "matrixUtils.h"

class threeStateModel : public replacementModel {
public:
	explicit threeStateModel(const MDOUBLE m1, const MDOUBLE m2, 
		const MDOUBLE m3, const MDOUBLE m4,const Vdouble &freq, bool useMarkovLimiting = true);
	threeStateModel(const threeStateModel& other) {*this = other;}
	virtual threeStateModel& operator=(const threeStateModel &other);
	virtual threeStateModel* clone() const { return new threeStateModel(*this); }
	virtual ~threeStateModel() {}
	const int alphabetSize() const {return 3;} // two states and an intermediate (both states at once)
	const MDOUBLE err_allow_for_pijt_function() const {return 1e-4;} // same as q2p definitions
	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const ;		
	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const{
		if (d==0.0)
			return _Q[i][j];
		errorMsg::reportError("Error in threeStateModel, dPij_dt called");
		return 0.0; // not supposed to be here
	}
	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const{
		errorMsg::reportError("Error in threeStateModel, d2Pij_dt2 called");
		return 0.0; // not supposed to be here
	}
	const MDOUBLE freq(const int i) const {
		if (i >= _freq.size()) 
			errorMsg::reportError("Error in threeStateModel::freq, i > size of frequency vector");
		return _freq[i];
	}
	const Vdouble getFreqs() const  {return _freq;}
	void setFreq(const Vdouble &freq);
	void setMu1(const MDOUBLE val) ;
	void setMu2(const MDOUBLE val) ;
	void setMu3(const MDOUBLE val) ;
	void setMu4(const MDOUBLE val) ;
	const MDOUBLE getMu1() const {return _gainState1;}
	const MDOUBLE getMu2() const {return _gainState0;}
	const MDOUBLE getMu3() const {return _lossState1;}
	const MDOUBLE getMu4() const {return _lossState0;}
	void computeMarkovLimitingDistribution(); // compute P(infinity), which specifies the stationary distribution

private:
	virtual void updateQ();
	void setEpsilonForZeroParams();
	bool checkIsNullModel();
	bool pijt_is_prob_value(MDOUBLE val) const;
	bool areFreqsValid(Vdouble freq) const; // tests if frequencies are valid (>0, sum=1)

private:
	
	MDOUBLE _gainState1; // _Q[0][2]
	MDOUBLE _gainState0; // _Q[1][2]
	MDOUBLE _lossState1; // _Q[2][0]
	MDOUBLE _lossState0; // _Q[2][1]
	VVdouble _Q;
	Vdouble _freq;
	bool _useMarkovLimiting; // should the markov limiting distribution be used to estimate the root frequencies
	mutable bool _bQchanged; //indicates whether the Q matrix was changed after the last Pij_t call
	mutable MDOUBLE _lastTcalculated;
	mutable VVdoubleRep _lastPtCalculated;


	
};

/*class gainLossModel : public replacementModel {
public:
explicit gainLossModel(const MDOUBLE m1, const MDOUBLE m2, const Vdouble freq);
virtual replacementModel* clone() const { return new gainLossModel(*this); }
gainLossModel(const gainLossModel& other): _q2pt(NULL) {*this = other;}
virtual gainLossModel& operator=(const gainLossModel &other);

virtual ~gainLossModel() {if (_q2pt) delete _q2pt;}
const int alphabetSize() const {return 3;} // two states and an intermediate (both states at once)
const MDOUBLE err_allow_for_pijt_function() const {return 1e-4;} // same as q2p definitions
const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const {
return _q2pt->Pij_t(i,j,d);
}
const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const{
return _q2pt->dPij_dt(i,j,d);
}
const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const{
return _q2pt->d2Pij_dt2(i,j,d);
}
const MDOUBLE freq(const int i) const {
if (i >= _freq.size()) 
errorMsg::reportError("Error in gainLossModel::freq, i > size of frequency vector");
return _freq[i];
}
void setMu1(const MDOUBLE val, bool isReversible=true) { _gainState1 = val; updateQ(isReversible);}
void setMu2(const MDOUBLE val,bool isReversible=true) { _gainState0 = val; updateQ(isReversible);}
const MDOUBLE getMu1() const {return _gainState1;}
const MDOUBLE getMu2() const {return _gainState0;}


protected:
virtual void updateQ(bool isReversible=true);
virtual void normalizeQ();


protected:
Vdouble _freq;
MDOUBLE _gainState1; 
MDOUBLE _gainState0; 
VVdouble _Q;
q2pt *_q2pt;



};
*/
/*
Q is a matrix of the following form:

0		1		01
0	1-m1	0		m1
1	0		1-m2	m2
01	(filled in assuming reversibility)

i.e. no direct change from state 0 to state 1 is allowed
*/

#endif // ___3STATE_MODEL


