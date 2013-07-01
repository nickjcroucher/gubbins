// $Id: khTest.cpp 962 2006-11-07 15:13:34Z privmane $

#include "definitions.h"
#include "logFile.h"
#include "errorMsg.h"
#include <cmath>

void makekhTest(const VVdouble & likelihoodVal, MDOUBLE diffNumOfFreeParam) {
	// assume that 2 trees are here.
	bool logValue = true;

	if (likelihoodVal.size() !=2) {errorMsg::reportError("errir un ");}

	const int n = likelihoodVal[0].size();
	
	MDOUBLE tmp1a = 0.0;
	MDOUBLE tmp1b = 0.0;
	MDOUBLE tmp1 = 0.0;
	MDOUBLE tmp2 = 0.0;
	MDOUBLE sum_k = 0.0;

	int k;
	for (k=0; k<n; ++k) {
		if (logValue==false) {
			tmp1a += log(likelihoodVal[1][k]);
			tmp1b += log(likelihoodVal[0][k]);
		}
		else {
			tmp1a += likelihoodVal[1][k];
			tmp1b += likelihoodVal[0][k];
		}
	}
	tmp1 = tmp1a-tmp1b;
	MDOUBLE difL = tmp1;
	for (k=0; k<n; ++k) {
		if (logValue==false) tmp2 = log(likelihoodVal[1][k])-log(likelihoodVal[0][k])-tmp1/static_cast<MDOUBLE>(n);
		else tmp2 = likelihoodVal[1][k]-likelihoodVal[0][k]-tmp1/static_cast<MDOUBLE>(n);
		sum_k += (tmp2*tmp2);
	}
	sum_k = sum_k * static_cast<MDOUBLE>(n) / static_cast<MDOUBLE>(n-1);
	LOG(1,<<" L1= "<<tmp1a<<" L2= "<<tmp1b<<endl);
	LOG(1,<<" delta L is "<<difL<<endl);
	LOG(1,<<" delta AIC is "<<difL-diffNumOfFreeParam<<endl);
	LOG(1,<<" var is "<<sum_k<<endl);
	LOG(1,<<" std is "<<sqrt(sum_k)<<endl);
	LOG(1,<<" z is "<<(difL -diffNumOfFreeParam )/sqrt(sum_k)<<endl);

// 	LOG(5,<<" L1= "<<tmp1a<<" L2= "<<tmp1b<<endl);
// 	LOG(5,<<" delta L is "<<difL<<endl);
// 	LOG(5,<<" delta AIC is "<<difL-diffNumOfFreeParam<<endl);
// 	LOG(5,<<" var is "<<sum_k<<endl;);
// 	LOG(5,<<" std is "<<sqrt(sum_k)<<endl);
// 	LOG(5,<<" z is "<<(difL -diffNumOfFreeParam )/sqrt(sum_k)<<endl);

}

