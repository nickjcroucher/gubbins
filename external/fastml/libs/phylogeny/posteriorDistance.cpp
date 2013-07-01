// $Id: posteriorDistance.cpp 5883 2009-02-06 10:42:11Z privmane $

#include "posteriorDistance.h"
#include "numRec.h"
#include "countTableComponent.h"
#include "likeDist.h"
#include "uniDistribution.h"
#include "someUtil.h"
#include "jcDistance.h"
#include <cmath>


class C_eval_gammaMLDistancesPosterior_d{ 
private:
    const stochasticProcess& _sp;
    const sequence& _s1;
    const sequence& _s2;
    const Vdouble* _weights;
    const VVdoubleRep& _posteriorProb; // pos, rate
public:
    C_eval_gammaMLDistancesPosterior_d(const stochasticProcess& sp,
									   const sequence& s1,
									   const sequence& s2,
									   const VVdoubleRep& posteriorProb,
									   const Vdouble  * weights)
		:  _sp(sp),
		   _s1(s1),
		   _s2(s2),
		   _weights(weights),
		   _posteriorProb(posteriorProb)
		{};


    MDOUBLE operator() (MDOUBLE dist) {
		MDOUBLE sumL=0.0;
		doubleRep posLikelihood = 0.0;
		MDOUBLE posLikelihood_d = 0.0;
		for (int pos=0; pos < _s1.seqLen(); ++pos){
			if (_s1.isUnknown(pos) && _s2.isUnknown(pos)) continue; // the case of two unknowns
			posLikelihood = 0.0;
			posLikelihood_d = 0.0;
			if (_s1.isUnknown(pos) && _s2.isSpecific(pos)) { 
				// this is the more complicated case, where s1 = ?, s2 = specific
				posLikelihood = _sp.freq(_s2[pos]);
				posLikelihood_d =0.0;
			}
			else if (_s2.isUnknown(pos) && _s1.isSpecific(pos)) {
				posLikelihood = _sp.freq(_s1[pos]);
				posLikelihood_d =0.0;
			} else {
				for (int rateCategor = 0; rateCategor<_sp.categories(); ++rateCategor) {
					MDOUBLE rate = _sp.rates(rateCategor);
					MDOUBLE pij= 0.0;
					MDOUBLE dpij=0.0;
					if (_s1.isSpecific(pos) && _s2.isSpecific(pos)) {//simple case, where AA i is changing to AA j
						pij= _sp.Pij_t(_s1[pos],_s2[pos],dist*rate);
						dpij= _sp.dPij_dt(_s1[pos],_s2[pos],dist*rate)*rate;
						doubleRep tmp =  _sp.freq(_s1[pos])*_posteriorProb[pos][rateCategor];
						posLikelihood += pij *tmp;
						posLikelihood_d += dpij*convert(tmp);
					} 
					else {// this is the most complicated case, when you have combinations of letters,
						// for example B in one sequence and ? in the other.
						for (int iS1 =0; iS1< _sp.alphabetSize(); ++iS1) {
							for (int iS2 =0; iS2< _sp.alphabetSize(); ++iS2) {
								if ((_s1.getAlphabet()->relations(_s1[pos],iS1)) &&
									(_s2.getAlphabet()->relations(_s2[pos],iS2))) {
									doubleRep exp = _sp.freq(iS1)*_posteriorProb[pos][rateCategor];;
									posLikelihood += exp* _sp.Pij_t(iS1,iS2,dist*rate);
									posLikelihood_d += convert(exp) * _sp.dPij_dt(iS1,iS2,dist*rate)*rate;
								}
							}
						}
					}
				}// end of for rate categories
			}
			assert(posLikelihood!=0.0);
			sumL += posLikelihood_d/convert(posLikelihood)*(_weights ? (*_weights)[pos]:1.0);
		}
		return -sumL;
    };
};

class C_eval_gammaMLDistancesPosterior{ 
private:
    const stochasticProcess& _sp;
    const sequence& _s1;
    const sequence& _s2;
    const Vdouble* _weights;
    const VVdoubleRep& _posteriorProb; // pos, rate
public:
    C_eval_gammaMLDistancesPosterior(const stochasticProcess& sp,
									 const sequence& s1,
									 const sequence& s2,
									 const VVdoubleRep& posteriorProb,
									 const Vdouble  * weights):  _sp(sp),
																 _s1(s1),
																 _s2(s2),
																 _weights(weights), 
																 _posteriorProb(posteriorProb)
		{};


    MDOUBLE operator() (MDOUBLE dist) {
		/*DEBUG LOG(9,<<"C_eval_gammaMLDistancesPosterior::operator():"); LOGDO(9,printTime(myLog::LogFile())); LOG(9,<<": dist = "<<dist<<endl); DEBUG*/
		MDOUBLE sumL=0.0;
		doubleRep posLikelihood = 0.0;

		for (int pos=0; pos < _s1.seqLen(); ++pos){
			/*DEBUG LOG(9,<<"C_eval_gammaMLDistancesPosterior::operator():"); LOGDO(9,printTime(myLog::LogFile())); LOG(9,<<": pos = "<<pos<<endl); DEBUG*/
			if (_s1.isUnknown(pos) && _s2.isUnknown(pos)) continue; // the case of two unknowns
						/*DEBUG LOG(9,<<"_posteriorProb  ="<<_posteriorProb[pos]<<endl); DEBUG*/
			posLikelihood = 0.0;
			/*DEBUG LOG(9,<<"posLikelihood = "<<posLikelihood<<endl); DEBUG*/
			if (_s1.isUnknown(pos) && _s2.isSpecific(pos)) { 
				// this is the more complicated case, where s1 = ?, s2 = specific
				posLikelihood = _sp.freq(_s2[pos]);
			}
			else if (_s2.isUnknown(pos) && _s1.isSpecific(pos)) {
				posLikelihood = _sp.freq(_s1[pos]);
			} else {
				for (int rateCategor = 0; rateCategor<_sp.categories(); ++rateCategor) {
					MDOUBLE rate = _sp.rates(rateCategor);
					/*DEBUG LOG(9,<<"rate = "<<rate<<endl); DEBUG*/
					MDOUBLE pij= 0.0;
					if (_s1.isSpecific(pos) && _s2.isSpecific(pos)) {//simple case, where AA i is changing to AA j
						/*DEBUG LOG(9,<<"Both are specific"<<endl); DEBUG*/
						pij= _sp.Pij_t(_s1[pos],_s2[pos],dist*rate);
						doubleRep exp =  _sp.freq(_s1[pos])*_posteriorProb[pos][rateCategor];
						/*DEBUG LOG(9,<<"exp = "<<exp<<endl); DEBUG*/
						posLikelihood += pij *exp;
						/*DEBUG LOG(9,<<"posLikelihood = "<<posLikelihood<<endl); DEBUG*/
					} 
					else {// this is the most complicated case, when you have combinations of letters,
						// for example B in one sequence and ? in the other.
						/*DEBUG LOG(9,<<"One or both are non-specific"<<endl); DEBUG*/
						for (int iS1 =0; iS1< _sp.alphabetSize(); ++iS1) {
							for (int iS2 =0; iS2< _sp.alphabetSize(); ++iS2) {
								if ((_s1.getAlphabet()->relations(_s1[pos],iS1)) &&
									(_s2.getAlphabet()->relations(_s2[pos],iS2))) {
									doubleRep exp = _sp.freq(iS1)*_posteriorProb[pos][rateCategor];
									posLikelihood += exp* _sp.Pij_t(iS1,iS2,dist*rate);
								}
							}
						}
						/*DEBUG LOG(9,<<"posLikelihood = "<<posLikelihood<<endl); DEBUG*/
					}
				}// end of for rate categories
			}
			assert(posLikelihood!=0.0);
			sumL += log(posLikelihood)*(_weights ? (*_weights)[pos]:1);
		}
		/*DEBUG LOG(9,<<"C_eval_gammaMLDistancesPosterior::operator():"); LOGDO(9,printTime(myLog::LogFile())); LOG(9,<<": returning "<<(-sumL)<<endl); DEBUG*/
		return -sumL;
    };
};

posteriorDistance::posteriorDistance(const stochasticProcess & sp,
									 const VVdoubleRep & posteriorProb,
									 const MDOUBLE toll,
									 const MDOUBLE maxPairwiseDistance) 
    :
    likeDist(sp,toll,maxPairwiseDistance),_posteriorProb(posteriorProb) 
{}

posteriorDistance::posteriorDistance(stochasticProcess & sp,
									 const VVdoubleRep & posteriorProb,
									 const MDOUBLE toll,
									 const MDOUBLE maxPairwiseDistance) 
    :
    likeDist(sp,toll,maxPairwiseDistance),_posteriorProb(posteriorProb) 
{}

posteriorDistance::posteriorDistance(const stochasticProcess & sp,
									 const MDOUBLE toll,
									 const MDOUBLE maxPairwiseDistance) 
    :
    likeDist(sp,toll,maxPairwiseDistance),_posteriorProb(0) 
{}


posteriorDistance::posteriorDistance(stochasticProcess & sp,
									 const MDOUBLE toll,
									 const MDOUBLE maxPairwiseDistance) 
    :
    likeDist(sp,toll,maxPairwiseDistance),_posteriorProb(0) 
{}

posteriorDistance::posteriorDistance(const posteriorDistance& other):
	likeDist(static_cast<likeDist>(other)), _posteriorProb(other._posteriorProb)
{}



// distance is computed based on the posterior probability
const MDOUBLE posteriorDistance::giveDistance(const sequence& s1,
											  const sequence& s2,
											  const Vdouble  * weights,
											  MDOUBLE* score) const 
{
	/*DEBUG LOG(9,<<"posteriorDistance::giveDistance - start"<<endl); LOGDO(9,printTime(myLog::LogFile())); DEBUG*/
	const MDOUBLE ax=0,	cx=_maxPairwiseDistance;
	MDOUBLE bx=_jcDist.giveDistance(s1,s2,weights,score)/*=1.0*/;
	if (!(bx==bx)) bx = 1.0;
	if (!(bx>0.0)) bx = 0.000001;
	MDOUBLE dist=-1.0;
	MDOUBLE resL = -dbrent(ax,bx,cx,
						   C_eval_gammaMLDistancesPosterior(_sp,s1,s2,_posteriorProb,weights),
						   C_eval_gammaMLDistancesPosterior_d(_sp,s1,s2,_posteriorProb,weights),
						   _toll,
						   &dist);
	if (score) *score = resL;
	return dist;
}

// =============================
// OBSOLETE: this function was moved to pairwiseGammaDistance.cpp
class C_evalAlphaForPairOfSeq{
private:
    const countTableComponentGam& _ctc;
    stochasticProcess& _sp;
    const MDOUBLE _branchL;
public:
    C_evalAlphaForPairOfSeq(const countTableComponentGam& ctc,
							const MDOUBLE branchL,
							stochasticProcess& sp):_ctc(ctc), _sp(sp), _branchL(branchL) {};

    MDOUBLE operator() (MDOUBLE alpha) {
		(static_cast<gammaDistribution*>(_sp.distr()))->setAlpha(alpha);
		C_evalLikeDist cev(_ctc,_sp);
		MDOUBLE L=cev(_branchL);
		LOG(10,<<"check alpha="<<alpha<<", bl="<<_branchL<<" gives "<<L<<endl);
		return L;
    };
};

// OBSOLETE: this function was moved to pairwiseGammaDistance.cpp
// returns the best alpha.
MDOUBLE optimizeAlphaFixedDist(stochasticProcess & sp,
							   const countTableComponentGam & ctc,
							   const MDOUBLE branchL,
							   const vector<MDOUBLE>  * weights,
							   MDOUBLE* score=NULL){ // changes sp.
    MDOUBLE bestA=0.0;
    MDOUBLE bestQ=0.0;
    const MDOUBLE upperBoundOnAlpha = 15.0;
    const MDOUBLE epsilonAlphaOptimization = 0.01;
    const MDOUBLE cx=upperBoundOnAlpha;// left, midle, right limit on alpha
    const MDOUBLE bx=cx*0.3;
    const MDOUBLE ax=0.0;

	
    bestQ = -brent(ax,bx,cx,
				   C_evalAlphaForPairOfSeq(ctc,branchL,sp),
				   epsilonAlphaOptimization,
				   &bestA);
    (static_cast<gammaDistribution*>(sp.distr()))->setAlpha(bestA);
    if (score) *score = bestQ;
    return bestA;
}



// OBSOLETE: this function was moved to pairwiseGammaDistance.cpp
class C_eval_gammaMLAlpha{ 
private:
    const stochasticProcess& _sp;
    const sequence& _s1;
    const sequence& _s2;
    const MDOUBLE _distance;
    const Vdouble* _weights;
    //  const VVdoubleRep& _posteriorProb; // pos, rate
public:
    C_eval_gammaMLAlpha(const stochasticProcess& sp,
						const sequence& s1,
						const sequence& s2,
						const MDOUBLE distance,
						//		      const VVdoubleRep& posteriorProb,
						const Vdouble  * weights):  _sp(sp),
													_s1(s1),
													_s2(s2),
													_distance(distance),
													_weights(weights) 
		//						  _posteriorProb(posteriorProb)
		{};

    // this cast is required as the distribution within the
    // stochasticProcess is kept as the parent "distribution" class that
    // knows nothing of Alpha
    void setAlpha(MDOUBLE alpha) {
		(static_cast<gammaDistribution*>(_sp.distr()))->setAlpha(alpha);
    }


    MDOUBLE operator() (MDOUBLE alpha) {
		setAlpha(alpha);
		MDOUBLE likelihood = likeDist::evalLikelihoodForDistance(_sp,_s1,_s2,_distance,_weights);
		LOG(11,<<"check alpha="<<alpha<<", bl="<<_distance<<" gives "<<likelihood<<endl);
		return -likelihood;
    };
} ;


// OBSOLETE: this function was moved to pairwiseGammaDistance.cpp
// returns the best alpha.
MDOUBLE optimizeAlphaFixedDist( const sequence& s1,
								const sequence& s2,
								stochasticProcess & sp,
								const MDOUBLE branchL,
								const vector<MDOUBLE>  * weights,
								MDOUBLE* score=NULL){ // changes sp.
    MDOUBLE bestA=0.0;
    MDOUBLE bestQ=0.0;
    const MDOUBLE upperBoundOnAlpha = 15.0;
    const MDOUBLE epsilonAlphaOptimization = 0.01;
    const MDOUBLE cx=upperBoundOnAlpha;// left, midle, right limit on alpha
    const MDOUBLE bx=cx*0.3;
    const MDOUBLE ax=0.0;

	
    bestQ = -brent(ax,bx,cx,
				   C_eval_gammaMLAlpha(sp,s1,s2,branchL,weights),
				   epsilonAlphaOptimization,
				   &bestA);
    (static_cast<gammaDistribution*>(sp.distr()))->setAlpha(bestA);
    if (score) *score = bestQ;
    return bestA;
}



MDOUBLE posteriorDistance::giveInitialGuessOfDistance(
    const sequence& s1,
    const sequence& s2,
    const vector<MDOUBLE>  * weights,
    MDOUBLE* score) const {
    uniDistribution ud;
    stochasticProcess uniSp(&ud,_sp.getPijAccelerator());
    likeDist ld(uniSp);
    return (ld.giveDistance(s1,s2,weights,score));
}

// OBSOLETE?  What's the difference between this function and giveDistanceOptAlphaForPairOfSequences???
MDOUBLE posteriorDistance::giveDistanceOptAlphaForEachPairOfSequences(	const sequence& s1,
																		const sequence& s2,
																		const vector<MDOUBLE>  * weights,
																		MDOUBLE* score,
																		MDOUBLE* alpha) const {

    MDOUBLE toll = 0.0001;
	
    MDOUBLE resL = 0.0;
    MDOUBLE resQ = 0.0;
    MDOUBLE currentDistance = giveInitialGuessOfDistance(s1,s2,weights,&resL);
	
    countTableComponentGam ctc; // from technical reasons.
    ctc.countTableComponentAllocatePlace(_sp.alphabetSize(),_sp.categories());
	
    stochasticProcess tmpSp(_sp);
    for (int z=0; z<s1.seqLen(); ++z) {
		for (int j=0; j < tmpSp.categories(); ++j) {
            ctc.addToCounts(s1[z],s2[z],j,weights?(*weights)[z]:tmpSp.ratesProb(j));
		}
    }
    const int maxIter = 30;
    MDOUBLE newDist = 0.0;
    MDOUBLE lastBestAlpha = 0.0;
    for (int i=0; i < maxIter; ++i) {
		lastBestAlpha = optimizeAlphaFixedDist(tmpSp,ctc,currentDistance,weights,&resL); // changes sp.
		(static_cast<gammaDistribution*>(tmpSp.distr()))->setAlpha(lastBestAlpha);
		LOG(8,<<"lastBestAlpha="<<lastBestAlpha<<"("<<(static_cast<gammaDistribution*>(tmpSp.distr()))->getAlpha()<<")"<<"\t L="<<resL<<"\t");
		likeDist tmpld(tmpSp); // we must create a new ld, that will include the stochastic process with the new alpha
		newDist = tmpld.giveDistance(ctc,resQ);
		LOG(8,<<"dist="<<newDist<<endl);
		if (fabs(newDist-currentDistance)<toll) break;
		currentDistance = newDist;
    }
    if (score) *score = resL;
    if (alpha) *alpha =   lastBestAlpha;
    assert (newDist >=0);
    return newDist;
	
}



// OBSOLETE: this function was moved to pairwiseGammaDistance.cpp
MDOUBLE posteriorDistance::giveDistanceOptAlphaForPairOfSequences(	const sequence& s1,
																	const sequence& s2,
																	const vector<MDOUBLE>  * weights,
																	MDOUBLE* score,
																	MDOUBLE* alpha) const {
  
    MDOUBLE toll = 0.0001;
	
    MDOUBLE resL = 0.0;
    MDOUBLE currentDistance = giveInitialGuessOfDistance(s1,s2,weights,&resL);
	
    countTableComponentGam ctc; // from technical reasons.
	
    stochasticProcess tmpSp(_sp);
	
    const int maxIter = 30;
    MDOUBLE newDist = 0.0;
    MDOUBLE lastBestAlpha = 0.0;
    for (int i=0; i < maxIter; ++i) {
		lastBestAlpha = optimizeAlphaFixedDist(s1, s2, tmpSp, currentDistance, weights, &resL); // changes sp.
		LOG(8,<<"lastBestAlpha="<<lastBestAlpha<<"("<<"\t L="<<resL<<"\t");
		likeDist tmpld(tmpSp); // we must create a new ld, that will include the stochastic process with the new alpha
		newDist = tmpld.giveDistance(s1, s2, weights, &resL);
		LOG(8,<<"dist="<<newDist<<"(L="<<resL<<")"<<endl);
		if (fabs(newDist-currentDistance)<toll) break;
		currentDistance = newDist;
    }
    if (score) *score = resL;
    if (alpha) *alpha = lastBestAlpha;
    assert (newDist >=0);
    return newDist;
	
}
