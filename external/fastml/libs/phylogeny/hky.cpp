// $Id: hky.cpp 4291 2008-06-23 10:23:10Z itaymay $

#include "hky.h"
#include "errorMsg.h"

hky::hky(const MDOUBLE inProb_a,
					const MDOUBLE inProb_c,
					const MDOUBLE inProb_g,
					const MDOUBLE inProb_t,
					const MDOUBLE TrTv)  {
	_freq.resize(4);
	_freq[0] = inProb_a;	_freq[1] = inProb_c;
	_freq[2] = inProb_g;	_freq[3] = inProb_t;
	initParams(TrTv);
}


hky::hky(vector<MDOUBLE> inProbs, const MDOUBLE TrTv) : _freq(inProbs)
{
	if (inProbs.size()!=4)
		errorMsg::reportError("hky::hky(vector<MDOUBLE> inProbs, const MDOUBLE TrTv) : the size of inProbs is not 4");
	initParams(TrTv);
}

void hky::initParams(MDOUBLE TrTv) // init _a, _b, _c, and _y by using _freq and TrTv
{
	MDOUBLE In_k = TrTv*2; // k is defined as alpha / beta.
	// In k2p Tr/Tv = alpha / 2*beta.

	_c = 2*(_freq[0]*_freq[2]+_freq[3]*_freq[1]);
	_y = 2*(_freq[0]+_freq[2])*(_freq[1]+_freq[3]);
	// c*_a + y*_b = 1;
	//_a/_b = k;
	_b = 1.0 / (_c*In_k+_y);
	_a = _b*In_k;
}

void hky::changeTrTv(const MDOUBLE TrTv){
	MDOUBLE In_k = TrTv*2; // k is defined as alpha / beta.
				  // In k2p Tr/Tv = alpha / 2*beta.
	_b = 1.0 / (_c*In_k+_y);
	_a = _b*In_k;
}

MDOUBLE hky::getTrTv() const {
	return (_a/(2.0*_b));
}

const MDOUBLE hky::Pij_t(const int i, const int j, const MDOUBLE t) const {
	const MDOUBLE &pa = _freq[0];
	const MDOUBLE &pc = _freq[1];
	const MDOUBLE &pg = _freq[2];
	const MDOUBLE &pt = _freq[3];
	const MDOUBLE py = pc+pt;
	const MDOUBLE pr = pa+pg;

	const MDOUBLE &b = _b;
	const MDOUBLE &a = _a;
	const MDOUBLE lamda3 = -(py*b+pr*a);
	const MDOUBLE lamda4 = -(py*a+pr*b);

	MDOUBLE term1=0.0;
	MDOUBLE term2=0.0;
	MDOUBLE term3=0.0;
	MDOUBLE termAll=0.0;
	switch (i) {
	case 0:
		switch (j) {
			case 0:
				term1 = pa;
				term2 = exp(-b*t)*(py)*pa/pr;
				term3 = pg*exp(t*lamda3)/pr;
				termAll = term1 + term2+term3;
				return termAll;
	
				break;
			case 1:
				termAll = pc - exp(-b*t)*pc;
				return termAll;
	
				break;
			case 2:
				term1 = pg;
				term2 = exp(-b*t)*py*pg/pr;
				term3 = -pg*exp(t*lamda3)/pr;
				termAll = term1 + term2+term3;
				return termAll;

				break;
			case 3:
				termAll = pt - exp(-b*t)*pt;
				return termAll;

				break;
		}
	break;
		
	case 1:
		switch (j) {
			case 0:
				termAll = pa - exp(-b*t)*pa;
				return termAll;
				break;
			case 1:
				term1 = pc;
				term2 = exp(-b*t)*pr*pc/py;
				term3 = pt*exp(t*lamda4)/py;
				termAll = term1 + term2+term3;
				return termAll;


				break;
			case 2:
				termAll = pg - exp(-b*t)*pg;
				return termAll;
				break;

			case 3:
				term1 = pt;
				term2 = exp(-b*t)*pr*pt/py;
				term3 = -pt*exp(t*lamda4)/py;
				termAll = term1 + term2 + term3;
				return termAll;

				break;
		}
	break;
				
	case 2:
		switch (j) {
			case 0:
				term1 = pa;
				term2 = exp(-b*t)*py*pa/pr;
				term3 = -pa*exp(t*lamda3)/pr;
				termAll = term1 + term2+term3;

				return termAll;
				break;
			case 1:
				termAll = pc - exp(-b*t)*pc;
				return termAll;
				break;
			case 2:
				term1 = pg;
				term2 = exp(-b*t)*py*pg/pr;
				term3 =  pa*exp(t*lamda3)/pr;
				termAll = term1 + term2 + term3;

				return termAll;
				break;

			case 3:
				termAll = pt - exp(-b*t)*pt;

				return termAll;
				break;
		}
	break;
	case 3:
		switch (j) {
			case 0:
				termAll = pa - exp(-b*t)*pa;
				return termAll;
				break;
			case 1:
				term1 = pc;
				term2 = exp(-b*t)*pr*pc/py;
				term3 = -pc*exp(t*lamda4)/py;
				termAll = term1 + term2+term3;
				return termAll;


				break;
			case 2:
				termAll = pg - exp(-b*t)*pg;
				return termAll;
				break;

			case 3:
				term1 = pt;
				term2 = exp(-b*t)*(pr)*pt/(py);
				term3 = pc*exp(t*lamda4)/(py);
				termAll = term1 + term2 + term3;
				return termAll;

				break;
		}
		break;

	}
	return -1;
}

const MDOUBLE hky::dPij_dt(const int i,const int j, const MDOUBLE t) const {
	const MDOUBLE &pa = _freq[0];
	const MDOUBLE &pc = _freq[1];
	const MDOUBLE &pg = _freq[2];
	const MDOUBLE &pt = _freq[3];
	const MDOUBLE py = pc+pt;
	const MDOUBLE pr = pa+pg;

	const MDOUBLE &b = _b;
	const MDOUBLE &a = _a;
	const MDOUBLE lamda3 = -(py*b+pr*a);
	const MDOUBLE lamda4 = -(py*a+pr*b);

	MDOUBLE term1, term2, term3,termAll;
	
	switch (i) {
	case 0:
		switch (j) {
			case 0://ok
				term1 = 0;
				term2 = exp(-b*t)*(py)*pa/pr;
				term2 *= -b;
				term3 = pg*exp(t*lamda3)/pr;
				term3*= lamda3;
				termAll = term1 + term2+term3;
				return termAll;
	
				break;
			case 1://ok
				termAll =  b* exp(-b*t)*pc;
				return termAll;
	
				break;
			case 2://ok
				term1 = 0;
				term2 = (-b)*exp(-b*t)*py*pg/pr;
				term3 = -pg*exp(t*lamda3)/pr;
				term3*=lamda3;
				termAll = term1 + term2+term3;
				return termAll;

				break;
			case 3://ok
				termAll = b*exp(-b*t)*pt;
				return termAll;

				break;
		}
	break;
		
	case 1:
		switch (j) {
			case 0://ok
				termAll = b*exp(-b*t)*pa;
				return termAll;
				break;
			case 1://ok
				term1 = 0;
				term2 = (-b)*exp(-b*t)*pr*pc/py;
				term3 = lamda4*pt*exp(t*lamda4)/py;
				termAll = term1 + term2+term3;
				return termAll;
				break;
			case 2://ok
				termAll = b*exp(-b*t)*pg;
				return termAll;
				break;
			case 3://ok
				term1 = 0;
				term2 = (-b)*exp(-b*t)*pr*pt/py;
				term3 = (lamda4)*(-pt)*exp(t*lamda4)/py;
				termAll = term1 + term2 + term3;
				return termAll;
				break;
		}
	break;
	case 2:
		switch (j) {
			case 0://ok
				term1 = 0;
				term2 = (-b)*exp(-b*t)*py*pa/pr;
				term3 = lamda3*(-pa)*exp(t*lamda3)/pr;
				termAll = term1 + term2+term3;
				return termAll;
				break;
			case 1://ok
				termAll = b*exp(-b*t)*pc;
				return termAll;
				break;
			case 2://ok
				term1 = 0;
				term2 = (-b)*exp(-b*t)*py*pg/pr;
				term3 =  lamda3*pa*exp(t*lamda3)/pr;
				termAll = term1 + term2 + term3;
				return termAll;
				break;
			case 3://ok
				termAll = b*exp(-b*t)*pt;
				return termAll;
				break;
		}
	break;
	case 3:
		switch (j) {
			case 0://ok
				termAll = b*exp(-b*t)*pa;
				return termAll;
				break;
			case 1://ok
				term1 = 0;
				term2 = (-b)*exp(-b*t)*pr*pc/py;
				term3 = lamda4*(-pc)*exp(t*lamda4)/py;
				termAll = term1 + term2+term3;
				return termAll;
				break;
			case 2://ok
				termAll = b* exp(-b*t)*pg;
				return termAll;
				break;
			case 3://ok
				term1 = 0;
				term2 = (-b)*exp(-b*t)*(pr)*pt/(py);
				term3 = (lamda4)*pc*exp(t*lamda4)/(py);
				termAll = term1 + term2 + term3;
				return termAll;
				break;
		}
		break;
	}
	return -1;
}

const MDOUBLE hky::d2Pij_dt2(const int i,const int j, const MDOUBLE t) const {
	const MDOUBLE &pa = _freq[0];
	const MDOUBLE &pc = _freq[1];
	const MDOUBLE &pg = _freq[2];
	const MDOUBLE &pt = _freq[3];
	const MDOUBLE py = pc+pt;
	const MDOUBLE pr = pa+pg;

	const MDOUBLE &b = _b;
	const MDOUBLE &a = _a;
	const MDOUBLE lamda3 = -(py*b+pr*a);
	const MDOUBLE lamda4 = -(py*a+pr*b);

	MDOUBLE term1, term2, term3,termAll;
	
	switch (i) {
	case 0:
		switch (j) {
			case 0://ok2
				term1 = 0;
				term2 = b*b*exp(-b*t)*(py)*pa/pr;
				term3 = lamda3*lamda3*pg*exp(t*lamda3)/pr;
				termAll = term1 + term2+term3;
				return termAll;
	
				break;
			case 1://ok2
				termAll =  -b*b* exp(-b*t)*pc;
				return termAll;
				break;
			case 2://ok2
				term1 = 0;
				term2 = b*b*exp(-b*t)*py*pg/pr;
				term3 = lamda3*lamda3*(-pg)*exp(t*lamda3)/pr;
				termAll = term1 + term2+term3;
				return termAll;
				break;
			case 3://ok2
				termAll = -b*b*exp(-b*t)*pt;
				return termAll;
				break;
		}
	break;
	case 1:
		switch (j) {
			case 0://ok2
				termAll = -b*b*exp(-b*t)*pa;
				return termAll;
				break;
			case 1://ok2
				term1 = 0;
				term2 = b*b*exp(-b*t)*pr*pc/py;
				term3 = lamda4*lamda4*pt*exp(t*lamda4)/py;
				termAll = term1 + term2+term3;
				return termAll;
				break;
			case 2://ok2
				termAll = -b*b*exp(-b*t)*pg;
				return termAll;
				break;
			case 3://ok2
				term1 = 0;
				term2 = b*b*exp(-b*t)*pr*pt/py;
				term3 = lamda4*lamda4*(-pt)*exp(t*lamda4)/py;
				termAll = term1 + term2 + term3;
				return termAll;
				break;
		}
	break;
	case 2:
		switch (j) {
			case 0://ok2
				term1 = 0;
				term2 = b*b*exp(-b*t)*py*pa/pr;
				term3 = lamda3*lamda3*(-pa)*exp(t*lamda3)/pr;
				termAll = term1 + term2+term3;
				return termAll;
				break;
			case 1://ok2
				termAll = -b*b*exp(-b*t)*pc;
				return termAll;
				break;
			case 2://ok2
				term1 = 0;
				term2 = b*b*exp(-b*t)*py*pg/pr;
				term3 =  lamda3*lamda3*pa*exp(t*lamda3)/pr;
				termAll = term1 + term2 + term3;
				return termAll;
				break;
			case 3://ok2
				termAll = -b*b*exp(-b*t)*pt;
				return termAll;
				break;
		}
	break;
	case 3:
		switch (j) {
			case 0://ok2
				termAll = -b*b*exp(-b*t)*pa;
				return termAll;
				break;
			case 1://ok2
				term1 = 0;
				term2 = b*b*exp(-b*t)*pr*pc/py;
				term3 = lamda4*lamda4*(-pc)*exp(t*lamda4)/py;
				termAll = term1 + term2+term3;
				return termAll;
				break;
			case 2://ok2
				termAll = -b*b* exp(-b*t)*pg;
				return termAll;
				break;
			case 3://ok2
				term1 = 0;
				term2 = b*b*exp(-b*t)*(pr)*pt/(py);
				term3 = lamda4*lamda4*pc*exp(t*lamda4)/(py);
				termAll = term1 + term2 + term3;
				return termAll;
				break;
		}
		break;
	}
	return -1;
}

const MDOUBLE hky::dPij_tdBeta(const int i, const int j, const MDOUBLE t) const {
	const MDOUBLE &pa = _freq[0];
	const MDOUBLE &pc = _freq[1];
	const MDOUBLE &pg = _freq[2];
	const MDOUBLE &pt = _freq[3];
	const MDOUBLE &py = pc+pt;
	const MDOUBLE &pr = pa+pg;

	const MDOUBLE &b = _b;
	const MDOUBLE &a = _a;
	const MDOUBLE &lamda3 = -(py*b+pr*a);
	const MDOUBLE &lamda4 = -(py*a+pr*b);

	MDOUBLE term2, term3,termAll;
	
	const MDOUBLE& dlamda3= -py+_y*pr/_c;
	const MDOUBLE& dlamda4= -pr+_y*py/_c;

	switch (i) {
		
	case 0:
		switch (j) {
			case 0:
				term2 = (-t)*exp(-b*t)*(py)*pa/pr;
				term3 = t*dlamda3*pg*exp(t*lamda3)/pr;
				termAll = term2+term3;
				return termAll;
	
				break;
			case 1:
				termAll = t* exp(-b*t)*pc;
				return termAll;
	
				break;
			case 2:
				term2 = (-t)*exp(-b*t)*py*pg/pr;
				term3 = t*dlamda3*(-pg)*exp(t*lamda3)/pr;
				termAll = term2+term3;
				return termAll;

				break;
			case 3:
				termAll = t* exp(-b*t)*pt;
				return termAll;

				break;
		}
	break;
		
	case 1:
		switch (j) {
			case 0:
				termAll = t* exp(-b*t)*pa;
				return termAll;
				break;
			case 1:
				term2 = (-t)*exp(-b*t)*pr*pc/py;
				term3 = t*dlamda4*pt*exp(t*lamda4)/py;
				termAll = term2+term3;
				return termAll;


				break;
			case 2:
				termAll = t* exp(-b*t)*pg;
				return termAll;
				break;

			case 3:
				term2 = (-t)*exp(-b*t)*pr*pt/py;
				term3 = t*dlamda4*(-pt)*exp(t*lamda4)/py;
				termAll = term2 + term3;
				return termAll;

				break;
		}
	break;
				
	case 2:
		switch (j) {
			case 0:
				term2 = (-t)*exp(-b*t)*py*pa/pr;
				term3 = t*dlamda3*(-pa)*exp(t*lamda3)/pr;
				termAll = term2+term3;

				return termAll;
				break;
			case 1:
				termAll = t*exp(-b*t)*pc;
				return termAll;
				break;
			case 2:
				term2 = (-t)*exp(-b*t)*py*pg/pr;
				term3 =  t*dlamda3*pa*exp(t*lamda3)/pr;
				termAll = term2 + term3;

				return termAll;
				break;

			case 3:
				termAll = t* exp(-b*t)*pt;
				return termAll;
				break;
		}
	break;
	case 3:
		switch (j) {
			case 0:
				termAll = t* exp(-b*t)*pa;
				return termAll;
				break;
			case 1:
				term2 = (-t)*exp(-b*t)*pr*pc/py;
				term3 = t*dlamda4*(-pc)*exp(t*lamda4)/py;
				termAll = term2+term3;
				return termAll;


				break;
			case 2:
				termAll = t* exp(-b*t)*pg;
				return termAll;
				break;

			case 3:
				term2 = (-t)*exp(-b*t)*(pr)*pt/(py);
				term3 = t*dlamda4*pc*exp(t*lamda4)/(py);
				termAll = term2 + term3;
				return termAll;

				break;
		}
		break;

	}
	return -1;
}

//Q[0][1] = freq[1]*_b	; Q[0][2] = freq[2]*_a	; Q[0][3] = freq[3]*_b;
//Q[1][0] = freq[0]*_b; 					; Q[1][2] = freq[2]*_b	; Q[1][3] = freq[3]*_a;
//Q[2][0] = freq[0]*_a;	Q[2][1] = freq[1]*_b	;	 			; Q[2][3] = freq[3]*_b;
//Q[3][0] = freq[0]*_b;	Q[3][1] = freq[1]*_a	; Q[3][2] = freq[2]*_b;

