// $Id: likelihoodComputationFactors.cpp 962 2006-11-07 15:13:34Z privmane $

#include "definitions.h"
#include "tree.h"
#include "computeUpAlg.h"
#include "likelihoodComputationFactors.h"
#include <cmath>
#include <cassert>

using namespace likelihoodComputation;

MDOUBLE likelihoodComputation::getLOG_LofPos(const int pos,
					  const tree& et,
					  const sequenceContainer& sc,
					  const stochasticProcess& sp,
					  const MDOUBLE gRate){ // when there is a global rate for this position
// using the pij of stochastic process rather than pre computed pij's...
	vector<MDOUBLE> factors;
	computeUpAlg cup;
	suffStatGlobalHomPos ssc;
	cup.fillComputeUpSpecificGlobalRateFactors(et,sc,pos,sp,ssc,gRate,factors);

	doubleRep tmp = 0.0;
	for (int let = 0; let < sp.alphabetSize(); ++let) {
		doubleRep tmpLcat=
				ssc.get(et.getRoot()->id(),let)*
				sp.freq(let);;
		assert(tmpLcat>=0);
		tmp+=tmpLcat;
	}
	return log(tmp)-factors[et.getRoot()->id()]*log(10.0);
}

