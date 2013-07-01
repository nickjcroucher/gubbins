// 	$Id: likeDist2USSRV.h 962 2006-11-07 15:13:34Z privmane $	
#ifndef ___LIKE_DIST_2_USSRV_H
#define ___LIKE_DIST_2_USSRV_H

#include "definitions.h"
#include "countTableComponent.h"
#include "distanceMethod.h"
#include "stochasticProcess.h"
#include "logFile.h"
#include "ussrvModel.h"
#include <cmath>
using namespace std;

class likeDist2USSRV : public distanceMethod {
public:
	explicit likeDist2USSRV(const ussrvModel& model,
					  const MDOUBLE toll =0.0001,
					  const MDOUBLE maxPairwiseDistance = 5.0) :  _model(model) ,_toll(toll),_maxPairwiseDistance(maxPairwiseDistance) 
	{}
  
	likeDist2USSRV (const likeDist2USSRV& other):  _model(other._model) ,_toll(other._toll),_maxPairwiseDistance(other._maxPairwiseDistance) {};
	virtual likeDist2USSRV* clone() const {return new likeDist2USSRV(*this);}

	// THIS FUNCTION DOES NOT RETURN THE LOG LIKELIHOOD IN RESQ, BUT RATHER "Q", THE CONTRIBUTION of this edge
    // TO THE EXPECTED LOG-LIKELIHOOD (SEE SEMPHY PAPER).
	// NEVERTHELESS, THE t that optimizes Q is the same t that optimizes log-likelihood.
	const MDOUBLE giveDistance(	const countTableComponentGam& ctcBase,
								const countTableComponentHom& ctcSSRV,
								MDOUBLE& resQ,
								const MDOUBLE initialGuess= 0.03) const; // initial guess
	
	
	// returns the estimated ML distance between the 2 sequences.
	// if score is given, it will be the log-likelihood.
	//!!!!!!!!!!!!!!TO DO @@@@
	const MDOUBLE giveDistance(const sequence& s1,
						const sequence& s2,
						const vector<MDOUBLE>  * weights,
						MDOUBLE* score=NULL) const { 
							LOG(4,<<"likeDist2USSRV:giveDistance : This method should never be used" << endl); 
							return 1;}

	const MDOUBLE giveDistanceBrent(const countTableComponentGam& ctcBase,
									const countTableComponentHom& ctcSSRV,	
									MDOUBLE& resL,
									MDOUBLE initialGuess) const;
	
private:
	const ussrvModel& _model;
	const MDOUBLE _toll;
	const MDOUBLE _maxPairwiseDistance;

};


class C_evalLikeDist2USSRV{
private:
	const countTableComponentGam& _ctcBase;
	const countTableComponentHom& _ctcSSRV;
	const ussrvModel& _model;
public:
	C_evalLikeDist2USSRV(const countTableComponentGam& ctcBase,
						 const countTableComponentHom& ctcSSRV,
						 const ussrvModel& model):_ctcBase(ctcBase),_ctcSSRV(ctcSSRV), _model(model) {};

	MDOUBLE operator() (MDOUBLE dist) {
		const MDOUBLE epsilonPIJ = 1e-10;
		MDOUBLE sumL=0.0;
		MDOUBLE pij;
		int categor, alph1,alph2;
		// base model
		const stochasticProcess& baseSp = _model.getBaseModel();
		
		for (alph1=0; alph1 < _ctcBase.alphabetSize(); ++alph1){
			for (alph2=0; alph2 <  _ctcBase.alphabetSize(); ++alph2){
				for (categor = 0; categor < baseSp.categories(); ++categor) {				
					MDOUBLE rate = baseSp.rates(categor);
					pij= baseSp.Pij_t(alph1,alph2,dist*rate);
					if (pij<epsilonPIJ) pij = epsilonPIJ;//SEE REMARK (1) FOR EXPLANATION
					sumL += _ctcBase.getCounts(alph1,alph2,categor)*(log(pij)-log(baseSp.freq(alph2)));//*_sp.ratesProb(rateCategor);// removed.
					
				}
			}
		}
		
		// ssrv model
		const stochasticProcessSSRV& ssrvSp = _model.getSSRVmodel();
		for (alph1=0; alph1 < _ctcSSRV.alphabetSize(); ++alph1){
			for (alph2=0; alph2 <  _ctcSSRV.alphabetSize(); ++alph2){
				pij = ssrvSp.Pij_t(alph1,alph2,dist);
				if (pij<epsilonPIJ) pij = epsilonPIJ;
				sumL+=_ctcSSRV.getCounts(alph1,alph2)*(log(pij)-log(ssrvSp.freq(alph2)));//*_sp.ratesProb(rateCategor);// removed.
			}
		}
		LOG(12,<<"check bl="<<dist<<" gives "<<sumL<<endl);

		return -sumL;
	}
};

// REMARK 1: THE LINE if if (pij<epsilonPIJ) pij = epsilonPIJ
// There are cases when i != j, and t!=0, and yet pij =0, because of numerical problems
// For these cases, it is easier to assume pij is very small, so that log-pij don't fly...

// @@@@ doesn't work
class C_evalLikeDist_d_2USSRV{ // derivative.
public:
  C_evalLikeDist_d_2USSRV(const countTableComponentGam& ctcBase,
						const countTableComponentHom& ctcSSRV,
						const ussrvModel& model)    : _ctcBase(ctcBase), _ctcSSRV(ctcSSRV),_model(model) {};

private:
	const  countTableComponentGam& _ctcBase;
	const  countTableComponentHom& _ctcSSRV;
	const ussrvModel& _model;

public:
	MDOUBLE operator() (MDOUBLE dist) {
		MDOUBLE	sumDL=0.0;
		MDOUBLE pij, dpij;
		int categor, alph1,alph2;
		// Base model
		const stochasticProcess& spBase = _model.getBaseModel();
		for (alph1=0; alph1 <  _ctcBase.alphabetSize(); ++alph1){
			for (alph2=0; alph2 <  _ctcBase.alphabetSize(); ++alph2){
				for (categor = 0; categor<_model.noOfCategor(); ++categor) {
					MDOUBLE rate = spBase.rates(categor);
					MDOUBLE pij= spBase.Pij_t(alph1,alph2,dist);
					MDOUBLE dpij= spBase.dPij_dt(alph1,alph2,dist);
					
					sumDL+= _ctcBase.getCounts(alph1,alph2,categor)*dpij
									*rate/pij; 
				}
			}
		}
		// SSRV model
		const stochasticProcessSSRV& spSSRV = _model.getSSRVmodel();
		for (alph1=0; alph1 <  _ctcSSRV.alphabetSize(); ++alph1){
			for (alph2=0; alph2 <  _ctcSSRV.alphabetSize(); ++alph2){
				pij= spSSRV.Pij_t(alph1,alph2,dist);
				dpij= spSSRV.dPij_dt(alph1,alph2,dist);
				sumDL+= _ctcSSRV.getCounts(alph1,alph2)*dpij/pij; //rate=1;
			}
		}

		LOG(8,<<"derivation = "<<-sumDL<<endl);
		return -sumDL;
	};
};

#endif // ___LIKE_DIST_2_USSRV_H

