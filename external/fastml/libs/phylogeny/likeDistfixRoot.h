// $Id: likeDistfixRoot.h 4470 2008-07-17 15:37:40Z cohenofi $

#ifndef ___LIKE_DIST_H_GL_FIX_ROOT
#define ___LIKE_DIST_H_GL_FIX_ROOT

#include "definitions.h"
#include "countTableComponent.h"
#include "distanceMethod.h"
#include "stochasticProcess.h"
#include "logFile.h"
#include "jcDistance.h"
#include "sequenceContainer.h"
#include "unObservableData.h"
#include <cmath>
using namespace std;

class likeDistfixRoot : public distanceMethod {
public:
    // WARNING: the stochasticProcess is NOT copied.  The same object is used
    explicit likeDistfixRoot(const stochasticProcess& sp,
		      const MDOUBLE toll =0.0001,
		      const MDOUBLE maxPairwiseDistance = 5.0,
			  unObservableData*  unObservableData_p=NULL)
	:  _sp(sp),_nonConstSpPtr(NULL),_toll(toll),_maxPairwiseDistance(maxPairwiseDistance),_unObservableData_p(unObservableData_p) {}

  likeDistfixRoot(const likeDistfixRoot& other)
	:  _sp(other._sp),_nonConstSpPtr(other._nonConstSpPtr),_toll(other._toll),_maxPairwiseDistance(other._maxPairwiseDistance),_jcDist(other._jcDist) {}

  virtual likeDistfixRoot* clone() const {return new likeDistfixRoot(*this);}
    // This constructor allows non-const stochasticProcess so that likeDistfixRoot will be able to change alpha, etc.
    explicit likeDistfixRoot(stochasticProcess& sp,
		      const MDOUBLE toll =0.0001,
		      const MDOUBLE maxPairwiseDistance = 5.0)
	:  _sp(sp),_nonConstSpPtr(&sp),_toll(toll),_maxPairwiseDistance(maxPairwiseDistance) {}

    // THIS FUNCTION DOES NOT RETURN THE LOG LIKELIHOOD IN RESQ, BUT RATHER "Q", THE CONTRIBUTION of this edge
    // TO THE EXPECTED LOG-LIKELIHOOD (SEE SEMPHY PAPER).
    // NEVERTHELESS, THE t that optimizes Q is the same t that optimizes log-likelihood.
    const MDOUBLE giveDistance(const vector<countTableComponentGam>& ctc,
			       MDOUBLE& resQ,
			       const MDOUBLE initialGuess= 0.03) const; // initial guess

    // given two sequences, it evaluates the log likelihood.
    MDOUBLE evalLogLikelihoodGivenDistance(const sequence& s1,
					   const sequence& s2,
					   const MDOUBLE dis2evaluate);

    // returns the estimated ML distance between the 2 sequences.
    // if score is given, it will be the log-likelihood.
    const MDOUBLE giveDistance(const sequence& s1,
			       const sequence& s2,
			       const vector<MDOUBLE>  * weights,
			       MDOUBLE* score=NULL) const;

    // this function creates a countTableComponent (ctc) from the two sequences.
    // it then computes the distance from this ctc.
    // THIS FUNCTION DOES NOT RETURN THE LOG LIKELIHOOD IN score, BUT RATHER "Q", THE CONTRIBUTION of this edge
    // TO THE EXPECTED LOG-LIKELIHOOD (SEE SEMPHY PAPER).
    // NEVERTHELESS, THE t that optimizes Q is the same t that optimizes log-likelihood.
    //MDOUBLE giveDistanceThroughCTC(const sequence& s1,
				//   const sequence& s2,
				//   const vector<MDOUBLE>  * weights,
				//   MDOUBLE* score=NULL) const;

    const MDOUBLE giveLikelihood(const sequence& s1,
				 const sequence& s2,
				 MDOUBLE distance,
				 const vector<MDOUBLE>  * weights=NULL) const;

    // return the stochasticProcess 
    const stochasticProcess& getStochasticProcess() const {return _sp;}
    stochasticProcess& getNonConstStochasticProcess();
    bool isTheInternalStochasticProcessConst() const {return !_nonConstSpPtr;}
    MDOUBLE getToll() const {return _toll;}
    MDOUBLE getMaxPairwiseDistance() const {return _maxPairwiseDistance;}

protected:
    const stochasticProcess &_sp;
    stochasticProcess *_nonConstSpPtr;
    const MDOUBLE _toll;
    const MDOUBLE _maxPairwiseDistance;
    jcDistance _jcDist;
	unObservableData*  _unObservableData_p;

private:
    const MDOUBLE giveDistanceBrent(	const vector<countTableComponentGam>& ctc,
					MDOUBLE& resL,
					const MDOUBLE initialGuess= 0.03) const; // initial guess
    const MDOUBLE giveDistanceNR(	const countTableComponentGam& ctc,
					MDOUBLE& resL,
					const MDOUBLE initialGuess= 0.03) const; // initial guess



public:
    static MDOUBLE evalLikelihoodForDistance(const stochasticProcess& sp,
						       const sequence& s1,
						       const sequence& s2,
						       const MDOUBLE dist,
						       const vector<MDOUBLE>  * weights=NULL);

};


class C_evallikeDistfixRoot{
private:
    const vector<countTableComponentGam>& _ctc;
    const stochasticProcess& _sp;
	unObservableData*  _unObservableData_p;
public:
    C_evallikeDistfixRoot(const vector<countTableComponentGam>& ctc,		// ctc[letterAtRoot][rate][alph][alph]
		   const stochasticProcess& inS1, unObservableData*  unObservableData_p=NULL)
		   :_ctc(ctc), _sp(inS1),_unObservableData_p(unObservableData_p) {};

    MDOUBLE operator() (MDOUBLE dist) 
	{
		//if(_plogLforMissingData){
		//	sequenceContainer scZero;
		//	gainLossAlphabet alph;
		//	scZero.startZeroSequenceContainerGL(_sc, alph);
		//	*_plogLforMissingData = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,scZero,*_sp);
		//}
		const MDOUBLE epsilonPIJ = 1e-10;
		MDOUBLE sumL=0.0;
		MDOUBLE sumLtemp=0.0;
		for (int letterAtRoot = 0; letterAtRoot < _sp.alphabetSize(); ++letterAtRoot){
			for (int alph1=0; alph1 < _ctc[0].alphabetSize(); ++alph1){
				for (int alph2=0; alph2 <  _ctc[0].alphabetSize(); ++alph2){
					for (int rateCategor = 0; rateCategor<_sp.categories(); ++rateCategor) {
						MDOUBLE rate = _sp.rates(rateCategor);
								
						MDOUBLE pij= _sp.Pij_t(alph1,alph2,dist*rate);
						if (pij<epsilonPIJ) pij = epsilonPIJ;//SEE REMARK (1) FOR EXPLANATION
						sumLtemp = _ctc[letterAtRoot].getCounts(alph1,alph2,rateCategor)*(log(pij)-log(_sp.freq(alph2)));//*_sp.ratesProb(rateCategor);// removed.
						//sumLtemp = _ctc[letterAtRoot].getCounts(alph1,alph2,rateCategor)*(log(pij)-log(_sp.freq(letterAtRoot)));//*_sp.ratesProb(rateCategor);// removed.
						//if(_unObservableData_p)
						//	sumL = sumL/(1- exp(_unObservableData_p->getlogLforMissingData()));	//???
						sumL += sumLtemp;
					}
				}
			}
		}
		if(_unObservableData_p)
			sumL = sumL/(1- exp(_unObservableData_p->getlogLforMissingData()));	//???
		LOG(12,<<"check bl="<<dist<<" gives "<<sumL<<endl);

		return -sumL;
    };
};

// REMARK 1: THE LINE if if (pij<epsilonPIJ) pij = epsilonPIJ
// There are cases when i != j, and t!=0, and yet pij =0, because of numerical problems
// For these cases, it is easier to assume pij is very small, so that log-pij don't fly...

class C_evalLikeDist_dGLfixRoot{ // derivative.
public:
    C_evalLikeDist_dGLfixRoot(const vector<countTableComponentGam>& ctc,
		     const stochasticProcess& inS1)    : _ctc(ctc), _sp(inS1) {};
private:
    const  vector<countTableComponentGam>& _ctc;
    const stochasticProcess& _sp;
public:
    MDOUBLE operator() (MDOUBLE dist) {
	MDOUBLE	sumDL=0.0;
	for (int letterAtRoot = 0; letterAtRoot < _sp.alphabetSize(); ++letterAtRoot){
		for (int alph1=0; alph1 <  _ctc[0].alphabetSize(); ++alph1){
			for (int alph2=0; alph2 <  _ctc[0].alphabetSize(); ++alph2){
				for (int rateCategor = 0; rateCategor<_sp.categories(); ++rateCategor) {
					MDOUBLE rate = _sp.rates(rateCategor);

					MDOUBLE pij= _sp.Pij_t(alph1,alph2,dist*rate);
					MDOUBLE dpij = _sp.dPij_dt(alph1,alph2,dist*rate);
					sumDL+= _ctc[letterAtRoot].getCounts(alph1,alph2,rateCategor)*dpij //*_sp.ratesProb(rateCategor) : removed CODE_RED
					*rate/pij;
				}
			}
		}//cerr<<"derivation = "<<-sumDL<<endl;
	}
	return -sumDL;
    };
};

class C_evalLikeDist_d2GLfixRoot{ // second derivative.
public:
    C_evalLikeDist_d2GLfixRoot(const countTableComponentGam& ctc,
		      const stochasticProcess& inS1)    : _ctc(ctc), _sp(inS1) {};
private:
    const  countTableComponentGam& _ctc;
    const stochasticProcess& _sp;
public:
    MDOUBLE operator() (MDOUBLE dist) {
	MDOUBLE	sumDL=0.0;
	for (int alph1=0; alph1 <  _ctc.alphabetSize(); ++alph1){
	    for (int alph2=0; alph2 <  _ctc.alphabetSize(); ++alph2){
		for (int rateCategor = 0; rateCategor<_sp.categories(); ++rateCategor) {
		    MDOUBLE rate = _sp.rates(rateCategor);

		    MDOUBLE pij= _sp.Pij_t(alph1,alph2,dist*rate);
		    MDOUBLE dpij = _sp.dPij_dt(alph1,alph2,dist*rate);
		    MDOUBLE d2pij = _sp.d2Pij_dt2(alph1,alph2,dist*rate);
		    sumDL+= rate*_ctc.getCounts(alph1,alph2,rateCategor)*
			(pij*d2pij - dpij *dpij )/(pij*pij);
		}
	    }
	}
	return -sumDL;
    };
};

#endif

