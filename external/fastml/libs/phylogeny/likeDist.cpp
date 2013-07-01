// $Id: likeDist.cpp 5956 2009-03-15 10:00:36Z adist $

#include "likeDist.h"
#include "numRec.h"
#include "someUtil.h"

stochasticProcess& likeDist::getNonConstStochasticProcess() {
  if (!_nonConstSpPtr) {
    errorMsg::reportError("likeDist::getNonConstStochasticProcess: Can't give non-const stochasticProcess because the stochasticProcess that was given to the constructor of this likeDist object was const");
  }
  return *_nonConstSpPtr;
}

// ======================= functors needed for the computations =============

class C_evalLikeDistDirect{
private:
    const stochasticProcess& _sp;
    const sequence& _s1;
    const sequence& _s2;
    const vector<MDOUBLE>  * _weights;
public:
    C_evalLikeDistDirect(const stochasticProcess& inS1,
			 const sequence& s1,
			 const sequence& s2,
			 const vector<MDOUBLE>  * weights): _sp(inS1),_s1(s1),_s2(s2),_weights(weights) {};

    MDOUBLE operator() (MDOUBLE dist) const {
	return -likeDist::evalLikelihoodForDistance(_sp,_s1,_s2,dist,_weights);
    }
};

MDOUBLE likeDist::evalLikelihoodForDistance(const stochasticProcess& sp,
					    const sequence& s1,
					    const sequence& s2,
					    const MDOUBLE dist,
					    const vector<MDOUBLE>  * weights)  {
    MDOUBLE sumL=0.0; // sum of log likelihoods
    MDOUBLE posLikelihood = 0.0; // likelihood of a specific position
    for (int pos=0; pos < s1.seqLen(); ++pos){
	if (s1.isUnknown(pos) && s2.isUnknown(pos)) continue; // the case of two unknowns
	posLikelihood = 0.0;
	if (s1.isUnknown(pos) && s2.isSpecific(pos)) { 
	    // this is the more complicated case, where s1 = ?, s2 = specific
	    posLikelihood = sp.freq(s2[pos]);
	} else if (s2.isUnknown(pos) && s1.isSpecific(pos)) {
	    posLikelihood = sp.freq(s1[pos]);
	} else {
	    for (int rateCategor = 0; rateCategor<sp.categories(); ++rateCategor) {
		MDOUBLE rate = sp.rates(rateCategor);
		MDOUBLE pij= 0.0;
		if (s1.isSpecific(pos) && s2.isSpecific(pos)) {//simple case, where AA i is changing to AA j
		    pij= sp.Pij_t(s1[pos],s2[pos],dist*rate);
		    posLikelihood += pij * sp.freq(s1[pos])*sp.ratesProb(rateCategor);
		} else {// this is the most complicated case, when you have
		    // combinations of letters, for example B in one
		    // sequence and ? in the other.
		    for (int iS1 =0; iS1< sp.alphabetSize(); ++iS1) {
			for (int iS2 =0; iS2< sp.alphabetSize(); ++iS2) {
			    if ((s1.getAlphabet()->relations(s1[pos],iS1)) &&
				(s2.getAlphabet()->relations(s2[pos],iS2))) {
				posLikelihood += sp.freq(iS1)*sp.Pij_t(iS1,iS2,dist*rate)*sp.ratesProb(rateCategor);
			    }
			}
		    }
		}
	    } // end of for on the rates
	}
	assert(posLikelihood!=0.0);
	sumL += log(posLikelihood)*(weights ? (*weights)[pos]:1.0);
    }
    return sumL;
};

class C_evalLikeDistDirect_d{ // derivative.
private:
    const stochasticProcess& _sp;
    const sequence& _s1;
    const sequence& _s2;
    const vector<MDOUBLE>  * _weights;
public:
    C_evalLikeDistDirect_d(const stochasticProcess& sp,
			   const sequence& s1,
			   const sequence& s2,
			   const vector<MDOUBLE>  * weights): _sp(sp),_s1(s1),_s2(s2),_weights(weights) {};

    MDOUBLE operator() (MDOUBLE dist) const {
	MDOUBLE sumL=0.0; // sum of log likelihoods
	MDOUBLE posLikelihood = 0.0; // likelihood of a specific position
	MDOUBLE posLikelihood_d = 0.0; // derivative of the likelihood at a specific position
	for (int pos=0; pos < _s1.seqLen(); ++pos){
	    if (_s1.isUnknown(pos) && _s2.isUnknown(pos)) continue; // the case of two unknowns
	    posLikelihood = 0.0;
	    posLikelihood_d = 0.0;
	    if (_s1.isUnknown(pos) && _s2.isSpecific(pos)) { 
		// this is the more complicated case, where s1 = ?, s2 = specific
		posLikelihood = _sp.freq(_s2[pos]);
		posLikelihood_d =0.0;
	    } else if (_s2.isUnknown(pos) && _s1.isSpecific(pos)) {
		posLikelihood = _sp.freq(_s1[pos]);
		posLikelihood_d =0.0;
	    } else {
		for (int rateCategor = 0; rateCategor<_sp.categories(); ++rateCategor) {
		    MDOUBLE rate = _sp.rates(rateCategor);
		    MDOUBLE pij= 0.0;
		    MDOUBLE dpij=0.0;
		    if (_s1.isSpecific(pos) && _s2.isSpecific(pos)) {
			//simple case, where AA i is changing to AA j
			pij= _sp.Pij_t(_s1[pos],_s2[pos],dist*rate);
			dpij= _sp.dPij_dt(_s1[pos],_s2[pos],dist*rate)*rate;
			MDOUBLE tmp =  _sp.freq(_s1[pos])*_sp.ratesProb(rateCategor);
			posLikelihood += pij *tmp;
			posLikelihood_d += dpij*tmp;
		    } else {// this is the most complicated case, when you have combinations of letters,
			// for example B in one sequence and ? in the other.
			for (int iS1 =0; iS1< _sp.alphabetSize(); ++iS1) {
			    for (int iS2 =0; iS2< _sp.alphabetSize(); ++iS2) {
				if ((_s1.getAlphabet()->relations(_s1[pos],iS1)) &&
				    (_s2.getAlphabet()->relations(_s2[pos],iS2))) {
				    MDOUBLE exp = _sp.freq(iS1)*_sp.ratesProb(rateCategor);
				    posLikelihood += exp* _sp.Pij_t(iS1,iS2,dist*rate);
				    posLikelihood_d += exp * _sp.dPij_dt(iS1,iS2,dist*rate)*rate;
				}
			    }
			}
		    }
		}// end of for rate categories
	    }
	    assert(posLikelihood>0.0);
	    sumL += (posLikelihood_d/posLikelihood)*(_weights ? (*_weights)[pos]:1.0);
	}
	return -sumL;
    };
};


// THIS FUNCTION EVALUATES THE LIKELIHOOD GIVEN THE DISTANCE
MDOUBLE likeDist::evalLogLikelihoodGivenDistance(const sequence& s1, const sequence& s2,
						 const MDOUBLE dis2evaluate) {
    C_evalLikeDistDirect Cev(_sp,s1,s2,NULL);
    return -Cev.operator ()(dis2evaluate);
}

MDOUBLE likeDist::giveDistanceThroughCTC(	const sequence& s1,
						const sequence& s2,
						const vector<MDOUBLE>  * weights,
						MDOUBLE* score) const {
    // only in the case of homogenous model - work through pairwise EM like
    countTableComponentGam ctc;
    if (_sp.categories() != 1) {
	errorMsg::reportError("this function only work for homogenous model.");
    }
    ctc.countTableComponentAllocatePlace(s1.getAlphabet()->size(),1);
    for (int i=0; i<s1.seqLen(); ++i) {
	ctc.addToCounts(s1[i],s2[i],0,weights?(*weights)[i]:1.0);
    }
    MDOUBLE resL =0;
    return giveDistance(ctc,resL);
}

const MDOUBLE likeDist::giveDistance(const countTableComponentGam& ctc,
				     MDOUBLE& resQ,
				     const MDOUBLE initialGuess) const {
    //return giveDistanceNR(ctc,resL,initialGuess);
    return giveDistanceBrent(ctc,resQ,initialGuess);
}

const MDOUBLE likeDist::giveDistanceBrent(const countTableComponentGam& ctc,
					  MDOUBLE& resL,
					  const MDOUBLE initialGuess) const {
    const MDOUBLE ax=0,bx=initialGuess,cx=_maxPairwiseDistance,tol=_toll;
    MDOUBLE dist=-1.0;
    resL = -dbrent(ax,bx,cx,
		   C_evalLikeDist(ctc,_sp,_unObservableData_p),
		   C_evalLikeDist_d(ctc,_sp),
		   tol,
		   &dist);
    return dist;
}

template <typename regF, typename dF>
MDOUBLE myNRmethod(MDOUBLE low, MDOUBLE current, MDOUBLE high, regF f,
		   dF df, const MDOUBLE tol, const int max_it, int & zeroFound) { // finding zero of a function.
    zeroFound = 1;
    MDOUBLE currentF = f(current);
    if (fabs(currentF)<tol) return current;
    MDOUBLE lowF = f(low);
    MDOUBLE highF = f(high);
    if (((lowF>0) && (highF>0)) || ((lowF<0) && (highF<0))) {// unable to find a zero
	zeroFound = 0;
	return 0;
    }
    if (lowF>0) {// fixing things to be in the right order.
	MDOUBLE tmp = low;
	low = high;
	high = tmp;
	tmp = lowF;
	lowF = highF;
	highF = tmp;
    }
    if (currentF>0) {
	high = current;
	highF = currentF;
    } else {
	low = current;
	lowF = currentF;
    } // now the zero is between current and either low or high.

    MDOUBLE currentIntervalSize = fabs(low-high);
    MDOUBLE oldIntervalSize = currentIntervalSize;

    // we have to decide if we do NR or devide the interval by two:
    // we want to check if the next NR step is within our interval
    // recall the the next NR guess is Xn+1 = Xn - f(Xn) / f(Xn+1)
    // So we want (current - currentF/currentDF) to be between low and high
    for (int i=0 ; i < max_it; ++i) {
	MDOUBLE currentDF = df(current);
	MDOUBLE newGuess = current - currentF/currentDF;
	if ((newGuess<low && newGuess> high) || (newGuess>low && newGuess< high)) {
	    // in this case we should do a NR step.
	    current = newGuess;
	    currentF = f(current);
	    if (currentF > 0){
		high = current;
		highF = currentF;
	    } else {
		low = current;
		lowF = currentF;
	    }

	    oldIntervalSize = currentIntervalSize;
	    currentIntervalSize =fabs (high-low);
	    if (currentIntervalSize < tol) {
		return current;
	    }
	    //LOG(5,<<"NR: low= "<<low<<" high= "<<high<<endl);
	}
	else { // bisection
	    oldIntervalSize = currentIntervalSize;
	    currentIntervalSize /= 2.0;
	    current = (low+high)/2.0;
	    currentF = f(current);
	    if (currentF > 0){
		high = current;
		highF = currentF;
	    } else {
		low = current;
		lowF = currentF;
	    }
	    //LOG(5,<<"BIS: low= "<<low<<" high= "<<high<<endl);
	    if (currentIntervalSize < tol) {
		return current;
	    }

	}
    }
    errorMsg::reportError("to many iterations in myNR function");
    return 0;
}

const MDOUBLE likeDist::giveDistanceNR(	const countTableComponentGam& ctc,
					MDOUBLE& resL,
					const MDOUBLE initialGuess) const {
    //change bx so that it will be the current branch length!
    const MDOUBLE ax=0,bx=initialGuess,cx=_maxPairwiseDistance,tol=_toll;
    //	LOG(5,<<"===================================================\n");
    MDOUBLE dist=-1.0;
    int zeroFound = 0;
    dist = myNRmethod(ax,bx,cx,
		      C_evalLikeDist_d(ctc,_sp),
		      C_evalLikeDist_d2(ctc,_sp),
		      tol,
		      100,
		      zeroFound);// max it for NR;
    if (zeroFound == 0) {// there was an error finding a zero
	dist = bx;
    }

    return dist;
}











/*




const MDOUBLE likeDist::giveDistance( // the NR version.
						const countTableComponentGam& ctc,
						MDOUBLE& resL) const {
	LOG(5,<<"=============="<<endl);
	MDOUBLE oldGuess=0.05; // move to parameters.
	if (oldGuess<0) oldGuess=0.05; // move up.
	int max_it = 100;
	MDOUBLE oldDist =0;
	MDOUBLE currentDist =oldGuess;
	MDOUBLE newDer =VERYBIG;
	MDOUBLE oldDer =VERYBIG;
	//const MDOUBLE ax=0,bx=1.0,cx=_maxPairwiseDistance,tol=_toll;
	for (int i=0; i < max_it; ++i){
		MDOUBLE	sumDL=0.0;
		MDOUBLE	sumDL2=0.0;
		for (int alph1=0; alph1 <  ctc.alphabetSize(); ++alph1){
			for (int alph2=0; alph2 <  ctc.alphabetSize(); ++alph2){
				for (int rateCategor = 0; rateCategor<_s1.categories(); ++rateCategor) {
					MDOUBLE rate = _s1.rates(rateCategor);

					MDOUBLE pij= _s1.Pij_t(alph1,alph2,currentDist*rate);
					MDOUBLE dpij = _s1.dPij_dt(alph1,alph2,currentDist*rate);
					MDOUBLE dpij2 = _s1.d2Pij_dt2(alph1,alph2,currentDist*rate);
					if (pij==0) {
						pij = 0.000000001;
						dpij = 0.000000001;
					}
					sumDL+= ctc.getCounts(alph1,alph2,rateCategor)*dpij 
									*rate/pij;
					sumDL2+= ctc.getCounts(alph1,alph2,rateCategor)*rate*(pij*dpij2-dpij *dpij)
									/(pij*pij);
				}
			}
		}
		oldDer = newDer;
		newDer = sumDL;
		LOG(5,<<"\ndistance = "<<currentDist<<endl);
		LOG(5,<<"derivation = "<<sumDL<<endl);
		LOG(5,<<"sec derivation = "<<sumDL2<<endl);
		oldDist = currentDist;
		if ((fabs(newDer) < fabs(oldDer)) && (sumDL2 < 0)) {
			currentDist = currentDist - newDer/sumDL2;
		}
		else {
			currentDist = currentDist / 2;
		}
		MDOUBLE epsilonForDeriv = 0.001;// move up
		if (fabs(newDer) < epsilonForDeriv) break;
		
	}

	return currentDist;
}*/

const MDOUBLE likeDist::giveDistance(const sequence& s1,
				     const sequence& s2,
				     const vector<MDOUBLE>  * weights,
				     MDOUBLE* score) const {
						 
    const MDOUBLE ax=0,	cx=_maxPairwiseDistance,tol=_toll;
    MDOUBLE bx=_jcDist.giveDistance(s1,s2,weights,score)/*=1.0*/;
    if (!(bx==bx)) bx = 1.0;  // safety check that the JC distance did not return nan (not a number)
    MDOUBLE dist=-1.0;
    MDOUBLE resL = -dbrent(ax,bx,cx,
			   C_evalLikeDistDirect(_sp,s1,s2,weights),
			   C_evalLikeDistDirect_d(_sp,s1,s2,weights),
			   tol,
			   &dist);
    if (score) *score = resL;
    return dist;
}

const MDOUBLE likeDist::giveLikelihood(const sequence& s1,
				       const sequence& s2,
				       MDOUBLE distance,
				       const vector<MDOUBLE>  * weights) const
{


    C_evalLikeDistDirect evalDis(_sp,s1,s2,weights);
    return -evalDis(distance);

}
