// $Id: gammaUtilities.cpp 962 2006-11-07 15:13:34Z privmane $

#include "gammaUtilities.h"
#include "logFile.h"
#include "errorMsg.h"
#include <cmath>


//gser: returns the incomplete Gamma function evaluated by its series representation
void gser(MDOUBLE *gamser, MDOUBLE a, MDOUBLE x, MDOUBLE *gln)
{
	//MDOUBLE gammln(MDOUBLE xx);
	
	int n;
	MDOUBLE sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) LOG(1,<<"x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		LOG(1,<<"Too many interations in routine gser");
		return;
	}
}

//gcf: returns the complement of the incomplete Gamma function evaluated by its continued fraction representation
void gcf(MDOUBLE *gammcf, MDOUBLE a, MDOUBLE x, MDOUBLE *gln)
{
	//MDOUBLE gammln(MDOUBLE xx);
	int i;
	MDOUBLE an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) LOG(1,<<"a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

//gammp(a, x): computes the incomplete Gamma function which is:
// 1/Gamma(a) * (the integral from 0 to x of (t^(a-1)*e^(-t)) dt)
//gammp can be computed in two different ways: by a series representation (gser(..)) 
//or by a continued fraction representation (gcf(..))
//gammp chooses to function will be used, according to the values of a and x 
MDOUBLE gammp(MDOUBLE a, MDOUBLE x)
{
	//void gcf(MDOUBLE *gammcf, MDOUBLE a, MDOUBLE x, MDOUBLE *gln);
	//void gser(MDOUBLE *gamser, MDOUBLE a, MDOUBLE x, MDOUBLE *gln);
	MDOUBLE gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) LOG(1,<<"Invalid arguments in routine gammp");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}



//I add////////////


MDOUBLE gammq(MDOUBLE a, MDOUBLE x)
{
	void gcf(MDOUBLE *gammcf, MDOUBLE a, MDOUBLE x, MDOUBLE *gln);
	void gser(MDOUBLE *gamser, MDOUBLE a, MDOUBLE x, MDOUBLE *gln);
	MDOUBLE gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) LOG(1,<<"Invalid arguments in routine gammp");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0 - gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}
/*************************************************************************
// this function computed the ln of the gamma function
// The Gamma funnction: Gamma(xx) = integral from 0 to infinity of (t^(xx-1)*e^(-t)) dt.
*************************************************************************/
MDOUBLE gammln(MDOUBLE xx)
{
	MDOUBLE x,y,tmp,ser;
	static MDOUBLE cof[6]={
		static_cast<MDOUBLE>(76.18009172947146),
		static_cast<MDOUBLE>(-86.50532032941677),
		static_cast<MDOUBLE>(24.01409824083091),
		static_cast<MDOUBLE>(-1.231739572450155),
		static_cast<MDOUBLE>(0.1208650973866179e-2),
		static_cast<MDOUBLE>(-0.5395239384953e-5)
	};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015f;
	for (j=0;j<6;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

//
MDOUBLE search_for_z_in_dis_with_any_beta(MDOUBLE alpha,MDOUBLE beta, MDOUBLE ahoson)
{
	return (search_for_z_in_dis_with_beta_1(alpha,ahoson)/beta);
}

MDOUBLE search_for_z_in_dis_with_beta_1(MDOUBLE alpha, MDOUBLE ahoson)
{
	if ( ahoson>1 || ahoson<0 ) errorMsg::reportError("Error in function search_for_z_in_dis_with_beta_1");
	MDOUBLE left=0;
	MDOUBLE right=99999.0;
	MDOUBLE tmp=5000.0;
	MDOUBLE results=0.0;

	for (int i=0;i<100000000 ; i++)
	{
		results=gammp(alpha,tmp);
		if (fabs(ahoson-results)<ERR_FOR_GAMMA_CALC) {
			return tmp;
		}
		if (results>ahoson) {
			right=tmp;
		}
		else left=tmp;
		tmp=(right+left)/2;
	}
	cout << "ERROR in search_for_z_in_dis_with_beta_1() Alpha is: "<< alpha <<endl;
	errorMsg::reportError("Error in function search_for_z_in_dis_with_beta_1 - first bonderi is 0");// also quit the program
	return 0;
}

MDOUBLE the_avarage_r_in_category_between_a_and_b(MDOUBLE left, MDOUBLE right, MDOUBLE alpha, MDOUBLE beta, int k)
{// and and b are the border of percentile k)
  MDOUBLE tmp;
  tmp= gammp(alpha+1,right*beta) - gammp(alpha+1,left*beta);
  tmp= (tmp*alpha/beta)*k;
  return tmp;
}
