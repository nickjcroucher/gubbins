// 	$Id: betaUtilities.cpp 962 2006-11-07 15:13:34Z privmane $	
#include "definitions.h"
#include "betaUtilities.h"
#include "gammaUtilities.h"
#include "logFile.h"
#include "errorMsg.h"
#include <cmath>

/******************************
	Computes the inverse of the beta CDF: given a prob. value, calculates the x for which 
	the integral over 0 to x of beta CDF = prob.
	Adapted from: 
	1. Majumder and Bhattacharjee (1973) App. Stat. 22(3) 411-414
	and the corrections:
	2. Cran et al. (1977) App. Stat. 26(1) 111-114
	3. Berry et al. (1990) App. Stat. 39(2) 309-310
	and another adaptation made in the code of Yang (tools.c)
****************************/
MDOUBLE inverseCDFBeta(MDOUBLE a, MDOUBLE b, MDOUBLE prob){
	if(a<0 || b<0 || prob<0 || prob>1)  {
		errorMsg::reportError("error in inverseCDFBeta,illegal parameter");
	}
	if (prob == 0 || prob == 1)
		return prob;
 
	int maxIter=100;
	MDOUBLE epsilonLow=1e-300;
	MDOUBLE fpu=3e-308;
            
	/****** changing the tail direction (prob=1-prob)*/
	bool tail=false;
	MDOUBLE probA=prob;
	if (prob > 0.5) {
		prob = 1.0 - prob; 
		tail = true;
		MDOUBLE tmp=a;
		a=b;
		b=tmp;
	}
	MDOUBLE lnBetaVal=betaln(a,b);
	MDOUBLE x; 
            
	/****** calculating chi square evaluator */        
	MDOUBLE r = sqrt(-log(prob * prob));
	MDOUBLE y = r - (2.30753+0.27061*r)/(1.+ (0.99229+0.04481*r) * r);
            
	MDOUBLE chiSquare = 1.0/(9.0 * b);
	chiSquare = b*2 * pow(1.0 - chiSquare + y * sqrt(chiSquare), 3.0);
//	MDOUBLE chiSquare2=gammq(b,prob/2.0); //chi square valued of prob with 2q df
	MDOUBLE T=(4.0*a+2.0*b-2)/chiSquare;
 
 
	/****** initializing x0 */
	if (a > 1.0 && b > 1.0) {
		r = (y * y - 3.) / 6.;
		MDOUBLE s = 1. / (a*2. - 1.);
		MDOUBLE t = 1. / (b*2. - 1.);
		MDOUBLE h = 2. / (s + t);
		MDOUBLE w = y * sqrt(h + r) / h - (t - s) * (r + 5./6. - 2./(3.*h));
		x = a / (a + b * exp(w + w));
	}
	else {
		if (chiSquare<0){
			x=exp((log(b*(1-prob))+lnBetaVal)/b);
		}
		else if (T<1){
			x=exp((log(prob*a)+lnBetaVal)/a);
		}
		else {
			x=(T-1.0)/(T+1.0);
		}
	}
            
	if(x<=fpu || x>=1-2.22e-16)  x=(prob+0.5)/2; // 0<x<1 but to avoid underflow a little smaller
 
	/****** iterating with a modified version of newton-raphson */
	MDOUBLE adj, newX=x, prev=0;
	MDOUBLE yprev = 0.;
	adj = 1.;
 
	MDOUBLE eps = pow(10., -13. - 2.5/(probA * probA) - 0.5/(probA *probA));
	eps = (eps>epsilonLow?eps:epsilonLow);
 
	for (int i=0; i<maxIter; i++) {
		y = incompleteBeta(a,b,x);
		y = (y - prob) *
			exp(lnBetaVal + (1.0-a) * log(x) + (1.0-b) * log(1.0 - x)); //the classical newton-raphson formula
		if (y * yprev <= 0) 
			prev = (fabs(adj)>fpu?fabs(adj):fpu);
		MDOUBLE g = 1;
		for (int j=0; j<maxIter; j++) {
			adj = g * y;
			if (fabs(adj) < prev) {
				newX = x - adj; // new x 
				if (newX >= 0. && newX <= 1.) {
					if (prev <= eps || fabs(y) <= eps)      return(tail?1.0-x:x);;
					if (newX != 0. && newX != 1.0)  break;
				}
			}
			g /= 3.;
		}
		if (fabs(newX-x)<fpu) 
			return (tail?1.0-x:x);;
		x = newX;
		yprev = y;
	}
	return (tail?1.0-x:x);
}


/******************************
	Computes the average r value in percentile k whose boundaries are leftBound and rightBound
****************************/
MDOUBLE computeAverage_r(MDOUBLE leftBound, MDOUBLE rightBound, MDOUBLE alpha, MDOUBLE beta, int k){
	MDOUBLE tmp;
	tmp= incompleteBeta(alpha+1,beta,rightBound) - incompleteBeta(alpha+1,beta,leftBound);
	tmp= (tmp*alpha/(alpha+beta))*k;
	return tmp;
}
/******************************
	Computes the integral from 0 to x over the beta CDF:
	(1/Beta(alpha,beta))x^(alpha-1)*(1-x)^(beta-1) where 
	Beta(a,b)=Gamma(a)*Gamma(b)/Gamma(a+b)
****************************/
MDOUBLE incompleteBeta(MDOUBLE alpha, MDOUBLE beta, MDOUBLE x){
	MDOUBLE tmp;
	if (x<0 || x>1) {
		LOG(5,<<"Error in function incompleteBeta : invalid x =  "<<x<<" alpha = "<<alpha<<" beta= "<<beta<<endl);
		errorMsg::reportError("Error in function incompleteBeta : invalid x");
	}
	if (x==0 || x==1) tmp=0.0;
	else tmp=exp(alpha*log(x)+beta*log(1-x)-betaln(alpha,beta));
	
	if (x<((alpha+1)/(alpha+beta+2))) return tmp*betacf(alpha,beta,x)/alpha;
	return 1-tmp*betacf(beta,alpha,1-x)/beta;
}
MDOUBLE betacf(MDOUBLE a, MDOUBLE b, MDOUBLE x){
	int m, m2;
	MDOUBLE aa,c,d,del,h,qab,qam,qap;
	qab = a+b;
	qap = a+1;
	qam = a-1;
	c=1;
	d=1-qab*x/qap;
	if (fabs(d)<FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for(m=1;m<=ITMAX;m++){
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d = 1.0+aa*d;
		if (fabs(d)<FPMIN) d = FPMIN;
		c=1.0 + aa/c;
		if (fabs(c)<FPMIN) c = FPMIN;
		d = 1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d = 1.0+aa*d;
		if (fabs(d)<FPMIN) d = FPMIN;
		c = 1.0 + aa/c;
		if (fabs(c)<FPMIN) c = FPMIN;
		d = 1.0/d;
		del = d*c;
		h*=del;
		if (fabs(del-1.0) <= EPS) break;
	}
	if (m > ITMAX) LOG(5,<<"Error in function betacf : alpha || beta big ||MAXIT small"<<endl);
	return h;
}

MDOUBLE betaln(MDOUBLE alpha, MDOUBLE beta){
	return gammln(alpha)+gammln(beta)-gammln(alpha+beta);
}

