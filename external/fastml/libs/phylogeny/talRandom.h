// $Id: talRandom.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___TAL_RANDOM
#define ___TAL_RANDOM

#include "definitions.h"
#include "logFile.h"
#include <cmath>
#include <cassert>
#include <ctime>

class RandintTal {
	unsigned long randx;
public:
	RandintTal(long s=0) {randx=s;}
	void seedTal(long s) {randx=s;}
	int absTal(int x) {return x&0x7fffffff;}
	static MDOUBLE maxTal() {return 2147483648.0;}
	int drawTal() {return randx = randx*1103515245+12345;}
	MDOUBLE fdrawTal() {return absTal(drawTal())/maxTal();} //random number between zero and 1
};

class talRandom {
public:
	// note the number you get is between 0 and entry not including entry!
	static MDOUBLE giveRandomNumberBetweenZeroAndEntry(MDOUBLE entry) {
		MDOUBLE tm=r.fdrawTal();
		return (tm * entry);
	}

	static bool flipCoin() {
		return ((talRandom::giveRandomNumberBetweenZeroAndEntry(1.0)-0.5)>0);
	}

	// note the number you get is between 0 and entry not including entry!
	static int giveIntRandomNumberBetweenZeroAndEntry(int entry)	{
		return (int)(giveRandomNumberBetweenZeroAndEntry(entry));
	}

	static void setSeed(const unsigned long seed) {
	  r.seedTal(seed);
	}

	static const MDOUBLE rand_gaussian(const MDOUBLE mean, const MDOUBLE variance) {
		const int N=100;
		static MDOUBLE X;
		X=0.0-N/2; /* set mean to 0 */
		for (int ri = 0;ri< N;ri++){
			//    X += 1.0*rand()/RAND_MAX;
			X += giveRandomNumberBetweenZeroAndEntry(1.0);
		}
	
		/* for uniform randoms in [0,1], mu = 0.5 and var = 1/12 */
		/* adjust X so mu = 0 and var = 1 */
		
		//  X = X * sqrt(12 / N);       /* adjust variance to 1 */
		//  cout <<X * sqrt(variance*12.0/N) + mean<<" ";
		MDOUBLE g = X * sqrt(variance*12.0/N) + mean;
		return (g);
	}

	static MDOUBLE  SampleGamma(MDOUBLE Alpha, MDOUBLE Beta) {
		MDOUBLE x= SampleGammaNorm(Alpha)/Beta;
		//LOG(700, << "SampleGamma(" << Alpha << " " << Beta << ") = " << x << "\n");
		return x;
	}
	static MDOUBLE  SampleGamma(MDOUBLE Alpha) {
	  MDOUBLE x= SampleGamma(Alpha, Alpha);
		//LOG(700, << "SampleGamma(" << Alpha << ") = " << x << "\n");
		return x;
	}
	static MDOUBLE rand_exp(const MDOUBLE mean) {
		return - mean * log(giveRandomNumberBetweenZeroAndEntry(1.0));//pg 64: Ross, Simulation 2nd.
	}

	static MDOUBLE giveRandomNumberBetweenTwoPoints(const MDOUBLE lower_point, const MDOUBLE upper_point) {
		MDOUBLE u = giveRandomNumberBetweenZeroAndEntry(upper_point - lower_point);
		return (u + lower_point);
	}


private:
	static RandintTal r;
	
	// Routine to generate a gamma random variable with unit scale (beta = 1)
	static MDOUBLE  SampleGammaNorm(MDOUBLE dblAlpha) {
		assert(dblAlpha > 0.0);
		if( dblAlpha < 1.0 ) return DblGammaLessThanOne(dblAlpha);
		else if( dblAlpha > 1.0 ) return DblGammaGreaterThanOne(dblAlpha);
		return -log(giveRandomNumberBetweenZeroAndEntry(1.0));
	}  
	static MDOUBLE DblGammaGreaterThanOne(MDOUBLE dblAlpha);
	static MDOUBLE  DblGammaLessThanOne(MDOUBLE dblAlpha);


};
#endif 

