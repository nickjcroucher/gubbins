// $Id: likeDist2Codon.h 4699 2008-08-14 14:19:46Z privmane $

#ifndef ___LIKE_DIST_2_CODON_H
#define ___LIKE_DIST_2_CODON_H

#include "definitions.h"
#include "countTableComponent.h"
#include "distanceMethod.h"
#include "stochasticProcess.h"
#include "logFile.h"
#include "wYangModel.h"
#include <cmath>
using namespace std;

class likeDist2Codon : public distanceMethod {
public:
	explicit likeDist2Codon(const vector<stochasticProcess>& spVec,
					  const MDOUBLE toll =0.0001,
					  const MDOUBLE maxPairwiseDistance = 2.0) :  _spVec(spVec) ,_toll(toll),_maxPairwiseDistance(maxPairwiseDistance) {
	}

	likeDist2Codon (const likeDist2Codon& other):  _spVec(other._spVec) ,_toll(other._toll),_maxPairwiseDistance(other._maxPairwiseDistance) {};
	virtual likeDist2Codon* clone() const {return new likeDist2Codon(*this);}

	// THIS FUNCTION DOES NOT RETURN THE LOG LIKELIHOOD IN RESQ, BUT RATHER "Q", THE CONTRIBUTION of this edge
    // TO THE EXPECTED LOG-LIKELIHOOD (SEE SEMPHY PAPER).
	// NEVERTHELESS, THE t that optimizes Q is the same t that optimizes log-likelihood.
	const MDOUBLE giveDistance(	const countTableComponentGam& ctc,
								MDOUBLE& resQ,
								const MDOUBLE initialGuess= 0.03) const; // initial guess
	
	
	// returns the estimated ML distance between the 2 sequences.
	// if score is given, it will be the log-likelihood.
	//!!!!!!!!!!!!!!TO DO
	const MDOUBLE giveDistance(const sequence& s1,
						const sequence& s2,
						const vector<MDOUBLE>  * weights,
						MDOUBLE* score=NULL) const { return 1;}

	const MDOUBLE giveDistanceBrent(	const countTableComponentGam& ctc,
										   MDOUBLE& resL,
										   const MDOUBLE initialGuess) const;
	
private:
	const vector<stochasticProcess>& _spVec;
	const MDOUBLE _toll;
	const MDOUBLE _maxPairwiseDistance;

};


class C_evalLikeDist2Codon{
private:
	const countTableComponentGam& _ctc;
	const vector<stochasticProcess>& _spVec;
public:
	C_evalLikeDist2Codon(const countTableComponentGam& ctc,
					const vector<stochasticProcess>& inS1):_ctc(ctc), _spVec(inS1) {};

	MDOUBLE operator() (MDOUBLE dist) {
		const MDOUBLE epsilonPIJ = 1e-10;
		MDOUBLE sumL=0.0;
		for (int alph1=0; alph1 < _ctc.alphabetSize(); ++alph1){
			for (int alph2=0; alph2 <  _ctc.alphabetSize(); ++alph2){
				for (int categor = 0; categor<_spVec.size(); ++categor) {				
					MDOUBLE pij= _spVec[categor].Pij_t(alph1,alph2,dist);
					if (pij<epsilonPIJ) pij = epsilonPIJ;//SEE REMARK (1) FOR EXPLANATION
					sumL += _ctc.getCounts(alph1,alph2,categor)*(log(pij)-log(_spVec[categor].freq(alph2)));//*_sp.ratesProb(rateCategor);// removed.
				}
			}
		}
	//	LOG(5,<<"check bl="<<dist<<" gives "<<sumL<<endl);

		return -sumL;
	};
};

// REMARK 1: THE LINE if if (pij<epsilonPIJ) pij = epsilonPIJ
// There are cases when i != j, and t!=0, and yet pij =0, because of numerical problems
// For these cases, it is easier to assume pij is very small, so that log-pij don't fly...

class C_evalLikeDist_d_2Codon{ // derivative.
public:
  C_evalLikeDist_d_2Codon(const countTableComponentGam& ctc,
				 const vector<stochasticProcess>& inS1)    : _ctc(ctc), _spVec(inS1) {};
private:
	const  countTableComponentGam& _ctc;
	const vector<stochasticProcess>& _spVec;
public:
	MDOUBLE operator() (MDOUBLE dist) {
		MDOUBLE	sumDL=0.0;
		for (int alph1=0; alph1 <  _ctc.alphabetSize(); ++alph1){
			for (int alph2=0; alph2 <  _ctc.alphabetSize(); ++alph2){
				for (int categor = 0; categor<_spVec.size(); ++categor) {
					MDOUBLE selection = static_cast<wYangModel*>(_spVec[categor].getPijAccelerator()->getReplacementModel())->getW();
					MDOUBLE pij= _spVec[categor].Pij_t(alph1,alph2,dist);
					MDOUBLE dpij = _spVec[categor].dPij_dt(alph1,alph2,dist);
					sumDL+= _ctc.getCounts(alph1,alph2,categor)*dpij //*_sp.ratesProb(rateCategor) : removed CODE_RED
									*selection/pij;
				}
			}
		}
		//LOG(5,<<"derivation = "<<-sumDL<<endl);
		return -sumDL;
	};
};

#endif

