// $Id: numRec.h 5790 2009-01-19 22:29:26Z rubi $

// version 1.00
// last modified 2 Nov 2002

#ifndef ___NUM_REC
#define ___NUM_REC

#include <cmath>
#include <cassert>
#include <iostream>
using namespace std;
#include "definitions.h"
#include "errorMsg.h"
#include "uniformDistribution.h"
#include "logFile.h"

//#define VERBOS
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

//========================== function brent =========================================
template <typename regF>
MDOUBLE brent(MDOUBLE ax, MDOUBLE bx, MDOUBLE cx, regF f, MDOUBLE tol,
	      MDOUBLE *xmin) {

    const int ITMAX  = 100;
    const MDOUBLE CGOLD = 0.3819660f;
    const MDOUBLE ZEPS = 1.0e-10f;

    int iter;
    MDOUBLE a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
    MDOUBLE e=0.0;

    a=(ax < cx ? ax : cx);
    b=(ax > cx ? ax : cx);
    x=w=v=bx;
    fw=fv=fx=f(x);
	LOG(10,<<"brent, f("<<x<<")="<<fx<<endl);
    for (iter=1;iter<=ITMAX;iter++) {
	xm=0.5*(a+b);
	tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
	if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
	    *xmin=x;
	    return fx;
	}
	if (fabs(e) > tol1) {
	    r=(x-w)*(fx-fv);
	    q=(x-v)*(fx-fw);
	    p=(x-v)*q-(x-w)*r;
	    q=2.0*(q-r);
	    if (q > 0.0) p = -p;
	    q=fabs(q);
	    etemp=e;
	    e=d;
	    if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
		d=CGOLD*(e=(x >= xm ? a-x : b-x));
	    else {
		d=p/q;
		u=x+d;
		if (u-a < tol2 || b-u < tol2)
		    d=SIGN(tol1,xm-x);
	    }
	} else {
	    d=CGOLD*(e=(x >= xm ? a-x : b-x));
	}
	u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
	fu=f(u);
	LOG(10,<<"brent, f("<<u<<")="<<fu<<endl);
	if (fu <= fx) {
	    if (u >= x) a=x; else b=x;
	    v=w;w=x;x=u;
	    fv=fw;fw=fx; fx=fu;
	} else {
	    if (u < x) a=u; else b=u;
	    if (fu <= fw || w == x) {
		v=w;
		w=u;
		fv=fw;
		fw=fu;
	    } else if (fu <= fv || v == x || v == w) {
		v=u;
		fv=fu;
	    }
	}
    }
    errorMsg::reportError(" too many iterations in function, brent. "); // also quit the program
    return -1;
}

// ===================================== function dbrent ========================================
/* The efficiency of this function for likelihood computations can be improved by replacing 
   functors regF f and dF df with one objects that preforms the likelihood computation once 
   and produces both L(t) and dL(t)/dt.  This object will provide methods:
   MDOUBLE f(MDOUBLE x)
   MDOUBLE df(MDOUBLE x)
*/	

#define ITMAX 100
#define ZEPS 1.0e-10
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);

template <typename regF, typename dF>
MDOUBLE dbrent(MDOUBLE ax, MDOUBLE bx, MDOUBLE cx, regF f,
	       dF df, MDOUBLE tol, MDOUBLE *xmin) {

    int iter,ok1,ok2;
    MDOUBLE a,b,d=0.0,d1,d2,du,dv,dw,dx,e=0.0;
    MDOUBLE fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

    a=(ax < cx ? ax : cx);
    b=(ax > cx ? ax : cx);
    //ensuring x is between a and b
    if (bx>b) { x=w=v=b;b=bx;}
    else if (bx<a) {x=w=v=a; a=bx;}
    else x=w=v=bx;
	
    fw=fv=fx=f(x);
    assert(fv==fv);// throw an exception if answer is nan = not a number.
    dw=dv=dx=df(x);

    for (iter=1;iter<=ITMAX;iter++) {
	xm=0.5*(a+b);
#ifdef VERBOS
	//if (iter>10) cout<<"iteration: "<<iter<<" xm = "<<xm<<" x= "<<x<<" a= "<<a<<" b= "<<b<<" fx= "<<fx<<endl;
#endif
	tol1=tol*fabs(x)+ZEPS;
	tol2=2.0*tol1;

	if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
	    *xmin=x;
	    return fx;
	}
	if (fabs(e) > tol1) { 
	    d1=2.0*(b-a);
	    d2=d1;
	    if (dw != dx) d1=(w-x)*dx/(dx-dw);
	    if (dv != dx) d2=(v-x)*dx/(dx-dv);
	    u1=x+d1;
	    u2=x+d2;
	    ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
	    ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
	    olde=e;
	    e=d;
	    if (ok1 || ok2) { 
		if (ok1 && ok2)
		    d=(fabs(d1) < fabs(d2) ? d1 : d2);
		else if (ok1)
		    d=d1;
		else
		    d=d2;
		if (fabs(d) <= fabs(0.5*olde)) {
		    u=x+d;
		    if (u-a < tol2 || b-u < tol2)
			d=SIGN(tol1,xm-x);
		} else {
		    d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
	    } else {
		d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	    }
	} else {
	    d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	}
	if (fabs(d) >= tol1) {
	    u=x+d;
	    fu=f(u);
	} else {
	    u=x+SIGN(tol1,d); 
	    if (u<ax) u=x; // MY LATEST ADDITION!
	    fu=f(u);
	    if (fu > fx) {
		*xmin=x;
		return fx;
	    }
	}
	du=df(u);
	if (fu <= fx) {
	    if (u >= x) a=x; else b=x;
	    MOV3(v,fv,dv, w,fw,dw)
		MOV3(w,fw,dw, x,fx,dx)
		MOV3(x,fx,dx, u,fu,du)
		} else {
		    if (u < x) a=u; else b=u; 
		    if (fu <= fw || w == x) {
			MOV3(v,fv,dv, w,fw,dw)
			    MOV3(w,fw,dw, u,fu,du)
			    } else if (fu < fv || v == x || v == w) {
				MOV3(v,fv,dv, u,fu,du)
				    }
		}

    }
    errorMsg::reportError("Too many iterations in routine dbrent"); // also quit the program
    return -1;
}

/*================================== function rtbis =========================================
//Using bisection, find the root of the function func known to lie between 
x1 and x2. The return value is the root will be refined until its accuracy is +- xacc 
*/
template <typename regF>
MDOUBLE rtbis(regF func,MDOUBLE x1, MDOUBLE x2, MDOUBLE xacc) {
    const int max_number_of_iter = 100;
	
    MDOUBLE f = func(x1);
    MDOUBLE fmid = func(x2);
    if (f*fmid >=0.0) {
	errorMsg::reportError(" error in function rtbis, root must be bracketed for bisection in rtbis ");
	// also quit the program
    }

    MDOUBLE dx, rtb;
    if (f<0.0) {
	dx = x2-x1;
	rtb = x1;
    }
    else {
	dx = x1-x2;
	rtb = x2;
    }


    for (int j=1; j <= max_number_of_iter; ++j) {
	dx *= 0.5;
	MDOUBLE xmid = rtb+dx; 
	MDOUBLE fmid = func(xmid);
	if (fmid <= 0.0) rtb = xmid;
	if ((fabs(dx) < xacc) || (fmid == 0.0)) return rtb;
    }
    errorMsg::reportError("Error in function rtbis..."); // also quit the program...
    return -1.0;
}

//Given a function func and an initial guessed range (x1,x2), the routine expands the range
//geometrically until a root is bracketed by the returned values x1 and x2 (in which case zbrac retruns true)
//or until the range becomes large unacceptably large (in which case zbrac return false).
template <typename regF>
bool zbrac(regF func, MDOUBLE &x1, MDOUBLE &x2) {
    const int NTRY=50;
    const MDOUBLE FACTOR= 1.6;
    int j;
    MDOUBLE f1,f2;

    if (x1 == x2) 
	errorMsg::reportError("Bad initial range in zbrac");
    f1 = func(x1);
    f2 = func(x2);
    for (j = 0; j < NTRY; j++) 
    {
	if (f1 * f2 < 0.0) 
	    return true;
	if (fabs(f1) < fabs(f2))
	    f1=func(x1 += FACTOR*(x1-x2));
	else
	    f2=func(x2 += FACTOR*(x2-x1));
    }
    return false;
}

// ================================ function brent new ======================================

int MyJacobi(VVdouble &Insym, VVdouble &RightEigenV, Vdouble &EigenValues);
MDOUBLE sign(MDOUBLE a,MDOUBLE b);
MDOUBLE pythag(const MDOUBLE a, const MDOUBLE b);
void houseHolder(VVdouble &mat,VVdouble &Q);
void tred2(VVdouble &a, Vdouble &d, Vdouble &e);
void QL(Vdouble &d, Vdouble &e, VVdouble &z);
void computeEigenSystem(VVdouble &symmetricMatrix,VVdouble &eigenVectros,Vdouble &diagonal);
MDOUBLE performKSTest(const uniformDistribution& empiricalDist, Vdouble& observedDist); // perform Kolomogorov-Smirnoff test
MDOUBLE computeProbForKS (const MDOUBLE QsParam); // function called only by performKSTest



#endif

