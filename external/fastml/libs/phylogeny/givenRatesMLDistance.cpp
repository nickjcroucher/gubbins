// 	$Id: givenRatesMLDistance.cpp 962 2006-11-07 15:13:34Z privmane $	
#include "givenRatesMLDistance.h"
#include "numRec.h"

class C_eval_likelihoodOfDistanceGivenRates{ 
private:
    const stochasticProcess& _sp;
    const sequence& _s1;
    const sequence& _s2;
    const Vdouble& _rates;
    const Vdouble* _weights;

public:
    C_eval_likelihoodOfDistanceGivenRates(const stochasticProcess& sp,
					  const sequence& s1,
					  const sequence& s2,
					  const Vdouble& rates,
					  const Vdouble  * weights)
	:  _sp(sp),_s1(s1),_s2(s2),_rates(rates),_weights(weights)
	{};

    MDOUBLE operator() (MDOUBLE dist) const {
	MDOUBLE sumL=0.0; // sum of log likelihoods
	MDOUBLE posLikelihood = 0.0; // likelihood of a specific position
	for (int pos=0; pos < _s1.seqLen(); ++pos){
	    if (_s1.isUnknown(pos) && _s2.isUnknown(pos)) continue; // the case of two unknowns
	    posLikelihood = 0.0;
	    if (_s1.isUnknown(pos) && _s2.isSpecific(pos)) { 
		// this is the more complicated case, where _s1 = ?, _s2 = specific
		posLikelihood = _sp.freq(_s2[pos]);
	    } else if (_s2.isUnknown(pos) && _s1.isSpecific(pos)) {
		posLikelihood = _sp.freq(_s1[pos]);
	    } else {
		MDOUBLE rate = _rates[pos];
		MDOUBLE pij= 0.0;
		if (_s1.isSpecific(pos) && _s2.isSpecific(pos)) {
		    // the simple case, where AA i is changing to AA j
		    pij= _sp.Pij_t(_s1[pos],_s2[pos],dist*rate);
		    posLikelihood += pij * _sp.freq(_s1[pos]);
		} else {// this is the most complicated case, when you have
		    // combinations of letters, for example B in one
		    // sequence and ? in the other.
		    for (int iS1 =0; iS1< _sp.alphabetSize(); ++iS1) {
			for (int iS2 =0; iS2< _sp.alphabetSize(); ++iS2) {
			    if ((_s1.getAlphabet()->relations(_s1[pos],iS1)) &&
				(_s2.getAlphabet()->relations(_s2[pos],iS2))) {
				posLikelihood += _sp.freq(iS1)*_sp.Pij_t(iS1,iS2,dist*rate);
			    }
			}
		    }
		}
	    }
	    assert(posLikelihood>0.0);
	    sumL += log(posLikelihood)*(_weights ? (*_weights)[pos]:1.0);
	}
	return -sumL;
    };
};

class C_eval_likelihoodOfDistanceGivenRates_d{ // derivative.
private:
    const stochasticProcess& _sp;
    const sequence& _s1;
    const sequence& _s2;
    const Vdouble& _rates;
    const Vdouble* _weights;

public:
    C_eval_likelihoodOfDistanceGivenRates_d(const stochasticProcess& sp,
					  const sequence& s1,
					  const sequence& s2,
					  const Vdouble& rates,
					  const Vdouble  * weights)
	:  _sp(sp),_s1(s1),_s2(s2),_rates(rates),_weights(weights)
	{};

    MDOUBLE operator() (MDOUBLE dist) const {
	MDOUBLE sumL=0.0; // sum of log likelihoods
	MDOUBLE posLikelihood = 0.0; // likelihood of a specific position
	MDOUBLE posLikelihood_d = 0.0; // derivative of the likelihood at a specific position
	for (int pos=0; pos < _s1.seqLen(); ++pos){
	    if (_s1.isUnknown(pos) && _s2.isUnknown(pos)) continue; // the case of two unknowns
	    posLikelihood = 0.0;
	    posLikelihood_d = 0.0;
	    if (_s1.isUnknown(pos) && _s2.isSpecific(pos)) { 
		// this is the more complicated case, where _s1 = ?, _s2 = specific
		posLikelihood = _sp.freq(_s2[pos]);
		posLikelihood_d =0.0;
	    } else if (_s2.isUnknown(pos) && _s1.isSpecific(pos)) {
		posLikelihood = _sp.freq(_s1[pos]);
		posLikelihood_d =0.0;
	    } else {
		MDOUBLE rate = _rates[pos];
		MDOUBLE pij= 0.0;
		MDOUBLE dpij=0.0;
		if (_s1.isSpecific(pos) && _s2.isSpecific(pos)) {
		    // the simple case, where AA i is changing to AA j
		    pij= _sp.Pij_t(_s1[pos],_s2[pos],dist*rate);
		    dpij= _sp.dPij_dt(_s1[pos],_s2[pos],dist*rate)*rate;
		    MDOUBLE tmp =  _sp.freq(_s1[pos]);
		    posLikelihood += pij *tmp;
		    posLikelihood_d += dpij*tmp;
		} else {// this is the most complicated case, when you have
		    // combinations of letters, for example B in one
		    // sequence and ? in the other.
		    for (int iS1 =0; iS1< _sp.alphabetSize(); ++iS1) {
			for (int iS2 =0; iS2< _sp.alphabetSize(); ++iS2) {
			    if ((_s1.getAlphabet()->relations(_s1[pos],iS1)) &&
				(_s2.getAlphabet()->relations(_s2[pos],iS2))) {
				MDOUBLE exp = _sp.freq(iS1);
				posLikelihood += exp* _sp.Pij_t(iS1,iS2,dist*rate);
				posLikelihood_d += exp * _sp.dPij_dt(iS1,iS2,dist*rate)*rate;
			    }
			}
		    }
		}
	    }
	    assert(posLikelihood>0.0);
	    sumL += (posLikelihood_d/posLikelihood)*(_weights ? (*_weights)[pos]:1.0);
	}
	return -sumL;
    };
};

const MDOUBLE givenRatesMLDistance::giveDistance(const sequence& s1,
						 const sequence& s2,
						 const vector<MDOUBLE>  * weights,
						 MDOUBLE* score) const
{
    const MDOUBLE ax=0,bx=1.0,cx=_maxPairwiseDistance;
    MDOUBLE dist=-1.0;
    MDOUBLE resL = -dbrent(ax,bx,cx,
			   C_eval_likelihoodOfDistanceGivenRates(_sp,s1,s2,_rates,weights),
			   C_eval_likelihoodOfDistanceGivenRates_d(_sp,s1,s2,_rates,weights),
			   _toll,
			   &dist);
    if (score) *score = resL;
    return dist;
};
