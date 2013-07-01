// 	$Id: normalDist.cpp 962 2006-11-07 15:13:34Z privmane $	
#include "normalDist.h"
#include <cmath>

/*
 This function evaluates the standard normal density function-N(0,1):
 integral from -infinity to x over exp(-.5t^2/sqrt(2pi)) (copied from the web) using 
  Milton Abramowiz and Irene A Stegun. 
  Handbook of Mathematical Functions. 
  National Bureau of Standards, 1964. 
 */
MDOUBLE Phi(MDOUBLE x)
{
	if (x>6.0) return 1;
	if (x<-6.0) return 0;
	MDOUBLE b1=0.31938153;
	MDOUBLE b2=-0.356563782;
	MDOUBLE b3=1.781477937;
	MDOUBLE b4=-1.821255978;
	MDOUBLE b5=1.330274429;
	MDOUBLE p=0.2316419;
	MDOUBLE c2=0.3989423;
	MDOUBLE a=fabs(x);
	MDOUBLE t=1.0/(1.0+a*p);
	MDOUBLE b=c2*exp((-x)*(x/2.0));
	MDOUBLE n=((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
	n=1.0-b*n;
	if (x<0.0) n=1.0-n;
	return n;
}

/*
 Computes the inverse normal distribution function (downloaded from the web)
 i.e. computes x when c=Phi(x)
 */
MDOUBLE normsinv(MDOUBLE p)
{
    if (p<EPSILON) return VERYSMALL;
    if ((1-p)<EPSILON)return VERYBIG;
    MDOUBLE x(0.0);
    MDOUBLE q, r;
    if ((0 < p )  && (p < P_LOW))
	{
        q = sqrt(-2*log(p));
        x = (((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
    }
    else
	{
        if ((P_LOW <= p) && (p <= P_HIGH))
		{
            q = p - 0.5;
            r = q*q;
            x = (((((A1*r+A2)*r+A3)*r+A4)*r+A5)*r+A6)*q /(((((B1*r+B2)*r+B3)*r+B4)*r+B5)*r+1);
        }
        else
		{
            if ((P_HIGH < p)&&(p < 1))
			{
				q = sqrt(-2*log(1-p));
				x = -(((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
			}
		}
    }
    return x;
}

 
