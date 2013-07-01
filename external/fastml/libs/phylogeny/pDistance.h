// $Id: pDistance.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___P_DISTANCE
#define ___P_DISTANCE

#include "definitions.h"
#include "distanceMethod.h"
/*********************************************************
p distance computes distance by counting number of differences and dividing by length of seq.
Weights are an input vector for giving additional weight to positions in the sequences.
*******************************************************/
class pDistance : public distanceMethod {
public:
	explicit pDistance(){}
	const MDOUBLE giveDistance(	const sequence& s1,
								const sequence& s2,
								const vector<MDOUBLE>  * weights,
								MDOUBLE* score=NULL) const {//score is not used here
		MDOUBLE p =0;
		if (weights == NULL) {
			for (int i = 0; i < s1.seqLen() ; ++i) if (s1[i] != s2[i]) p++;
			p = p/s1.seqLen();
		} else {
			MDOUBLE len=0;
			for (int i = 0; i < s1.seqLen() ; ++i) {
				len +=((*weights)[i]);
				if (s1[i] != s2[i]) p+=((*weights)[i]);
			}
			p = p/len;
		}
		return p;
	}
  virtual   pDistance* clone() const {return new pDistance(*this);}

};

#endif
