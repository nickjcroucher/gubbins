// $Id: talRandom.cpp 962 2006-11-07 15:13:34Z privmane $

#include "talRandom.h"

RandintTal talRandom::r = static_cast<long>(time(0)) ;

MDOUBLE  talRandom::DblGammaGreaterThanOne(MDOUBLE dblAlpha) {
  // Code adopted from David Heckerman
  //-----------------------------------------------------------
  //	DblGammaGreaterThanOne(dblAlpha)
  //
  //	routine to generate a gamma random variable with unit scale and
  //      alpha > 1
  //	reference: Ripley, Stochastic Simulation, p.90 
  //	Chang and Feast, Appl.Stat. (28) p.290
  //-----------------------------------------------------------
    MDOUBLE rgdbl[6];
    
    rgdbl[1] = dblAlpha - 1.0;
    rgdbl[2] = (dblAlpha - (1.0 / (6.0 * dblAlpha))) / rgdbl[1];
    rgdbl[3] = 2.0 / rgdbl[1];
    rgdbl[4] = rgdbl[3] + 2.0;
    rgdbl[5] = 1.0 / sqrt(dblAlpha);
    
    for (;;)
      {
	MDOUBLE  dblRand1;
	MDOUBLE  dblRand2;
	do
	  {
	    dblRand1 = giveRandomNumberBetweenZeroAndEntry(1.0);
	    dblRand2 = giveRandomNumberBetweenZeroAndEntry(1.0);
	    
	    if (dblAlpha > 2.5)
	      dblRand1 = dblRand2 + rgdbl[5] * (1.0 - 1.86 * dblRand1);
	    
	  } while (!(0.0 < dblRand1 && dblRand1 < 1.0));
	
	MDOUBLE dblTemp = rgdbl[2] * dblRand2 / dblRand1;
	
	if (rgdbl[3] * dblRand1 + dblTemp + 1.0 / dblTemp <= rgdbl[4] ||
	    rgdbl[3] * log(dblRand1) + dblTemp - log(dblTemp) < 1.0)
	  {
	    return dblTemp * rgdbl[1];
	  }
      }
    assert(false);
    return 0.0;
}

MDOUBLE  talRandom::DblGammaLessThanOne(MDOUBLE dblAlpha){
//routine to generate a gamma random variable with 
//unit scale and alpha < 1
//reference: Ripley, Stochastic Simulation, p.88 
	MDOUBLE dblTemp;
	const MDOUBLE	dblexp = exp(1.0);
	for (;;){
		MDOUBLE dblRand0 = giveRandomNumberBetweenZeroAndEntry(1.0);
		MDOUBLE dblRand1 = giveRandomNumberBetweenZeroAndEntry(1.0);
		if (dblRand0 <= (dblexp / (dblAlpha + dblexp))){
			dblTemp = pow(((dblAlpha + dblexp) * dblRand0) /
			dblexp, 1.0 / dblAlpha);
			if (dblRand1 <= exp(-1.0 * dblTemp)) return dblTemp;
		} else {
			dblTemp = -1.0 * log((dblAlpha + dblexp) * (1.0 - dblRand0) /
			 (dblAlpha * dblexp)); 
			if (dblRand1 <= pow(dblTemp,dblAlpha - 1.0)) return dblTemp;
		}
	}
    assert(false);
    return 0.0;
}  // DblGammaLessThanOne
  
