// $Id: jcDistance.h 1928 2007-04-04 16:46:12Z privmane $

#ifndef ___JC_DISTANCE
#define ___JC_DISTANCE

#include "definitions.h"
#include "distanceMethod.h"
#include <typeinfo>
#include <cmath>
/*********************************************************
Jukes-Cantor distance method. 
Assumes no constraints on replacement from one state to another.
Receives size of alphabet in constructor, and this enables 
to have one class for JC-distance for nucleotides, a.a., and codons  
Weights are an input vector for giving additional weight to positions in the sequences.
*******************************************************/
class jcDistance : public distanceMethod {

public:
	explicit jcDistance() {}
	virtual jcDistance* clone() const{ return new jcDistance(*this);}

	const MDOUBLE giveDistance(	const sequence& s1,
								const sequence& s2,
								const vector<MDOUBLE>  * weights,
								MDOUBLE* score=NULL) const {//score is not used here

		if (typeid(s1.getAlphabet()) != typeid(s2.getAlphabet()))
			errorMsg::reportError("Error in jcDistance::giveDistance, s1 and s2 contain different type of alphabet");
		
		// pS1Base and pS2Base are references to s1 and s2 respectively. 
		// The method uses seq1 and seq2 and not s1 and s2, because when
		// the sequences contain mulAlphabet we must first convert them to the base alphabet
		const sequence* pS1Base(&s1);
		const sequence* pS2Base(&s2);
		const alphabet* alph = s1.getAlphabet();
		// if s1 and contains mulAlphabet
		const mulAlphabet* mulAlph = dynamic_cast<const mulAlphabet*>(alph);
		if (mulAlph!=NULL) {
			pS1Base = new sequence(s1,mulAlph->getBaseAlphabet());
			pS2Base = new sequence(s2,mulAlph->getBaseAlphabet());
		}
		
		int alphabetSize = pS1Base->getAlphabet()->size();
		
		//		const MDOUBLE MAXDISTANCE=2.0;
		const MDOUBLE MAXDISTANCE=15;
		
		MDOUBLE p =0;
		MDOUBLE len=0.0;
		if (weights == NULL) {
			for (int i = 0; i < pS1Base->seqLen() ; ++i) {
				if ((*pS1Base)[i]<0 || (*pS2Base)[i]<0) continue; //gaps and missing data.
				len+=1.0;
				if ((*pS1Base)[i] != (*pS2Base)[i]) p++;
			}
			if (len==0) p=1;
			else p = p/len;
		} else {
			for (int i = 0; i < pS1Base->seqLen() ; ++i) {
				if ((*pS1Base)[i]<0 || (*pS2Base)[i]<0) continue; //gaps and missing data.
				len += (*weights)[i];
				if ((*pS1Base)[i] != (*pS2Base)[i])  p+=((*weights)[i]);
			}
			if (len==0) p=1;
			else {
				p = p/len;
			}
		}
		if (pS1Base != &s1) {
			delete pS1Base;
			delete pS2Base;
		}

		const MDOUBLE inLog = 1 - (MDOUBLE)alphabetSize*p/(alphabetSize-1.0);
		if (inLog<=0) {
//			LOG(6,<<" DISTANCES FOR JC DISTANCE ARE TOO BIG");
//			LOG(6,<<" p="<<p<<endl);
			return MAXDISTANCE;
		}
		MDOUBLE dis = -1.0 * (1.0 - 1.0/alphabetSize) * log (inLog);
		return dis;
	}
};

class jcDistanceOLD : public distanceMethod {
// in this version, if you have
// a gap in front of a letter - it will be taken as a different
// and also the length of the pairwise comparison will be increased.
// in case of a gap-gap, it won't be a difference, but the length will
// be increase.

private:
	const int _alphabetSize;

public:
	explicit jcDistanceOLD(const int alphabetSize) : _alphabetSize(alphabetSize) {
	}
	explicit jcDistanceOLD(const jcDistanceOLD& other) : _alphabetSize(other._alphabetSize) {
	}
	virtual jcDistanceOLD* clone() const{ return new jcDistanceOLD(*this);}

	const MDOUBLE giveDistance(	const sequence& s1,
								const sequence& s2,
								const vector<MDOUBLE>  * weights,
								MDOUBLE* score=NULL) const {//score is not used here
//		const MDOUBLE MAXDISTANCE=2.0;
		const MDOUBLE MAXDISTANCE=15;
		
		MDOUBLE p =0;
		MDOUBLE len=0.0;
		if (weights == NULL) {
			for (int i = 0; i < s1.seqLen() ; ++i) {
				//if (s1[i]<0 || s2[i]<0) continue; //gaps and missing data.
				len+=1.0;
				if (s1[i] != s2[i]) p++;
			}
			if (len==0) p=1;
			else p = p/len;
		} else {
			for (int i = 0; i < s1.seqLen() ; ++i) {
				//if (s1[i]<0 || s2[i]<0) continue; //gaps and missing data.
				len += (*weights)[i];
				if (s1[i] != s2[i])  p+=((*weights)[i]);
			}
			if (len==0) p=1;
			else {
				p = p/len;
			}
		}
		const MDOUBLE inLog = 1 - (MDOUBLE)_alphabetSize*p/(_alphabetSize-1.0);
		if (inLog<=0) {
//			LOG(6,<<" DISTANCES FOR JC DISTANCE ARE TOO BIG");
//			LOG(6,<<" p="<<p<<endl);
			return MAXDISTANCE;
		}
		MDOUBLE dis = -1.0 * (1.0 - 1.0/_alphabetSize) * log (inLog);
		return dis;
	}
};
#endif
