// $Id: GLaguer.cpp 962 2006-11-07 15:13:34Z privmane $
#include "definitions.h"
#include "GLaguer.h"

#include "errorMsg.h"
#include "gammaUtilities.h"



GLaguer::GLaguer(const int pointsNum, const MDOUBLE alf, Vdouble & points, Vdouble & weights)
{
	gaulag(_points, _weights, alf, pointsNum);

	weights = _weights;
	points = _points;
}


//Input: alf = the alpha parameter of the Laguerre polynomials
//		 pointsNum = the polynom order
//Output: the abscissas and weights are stored in the vecotrs x and w, respectively. 
//Discreption: given alf, the alpha parameter of the Laguerre polynomials, the function returns the abscissas and weights
//			   of the n-point Guass-Laguerre quadrature formula.
//			   The smallest abscissa is stored in x[0], the largest in x[pointsNum - 1].
void GLaguer::gaulag(Vdouble &x, Vdouble  &w, const MDOUBLE alf, const int pointsNum)
{
	x.resize(pointsNum, 0.0);
	w.resize(pointsNum, 0.0);
	const int MAXIT=10000;
	const MDOUBLE EPS=1.0e-6;
	int i,its,j;
	MDOUBLE ai,p1,p2,p3,pp,z=0.0,z1;

	int n= x.size();
	for (i=0;i<n;i++) {
		//loops over the desired roots
		if (i == 0) { //initial guess for the smallest root
			z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		} else if (i == 1) {//initial guess for the second smallest root
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		} else { //initial guess for the other roots
			ai=i-1;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
				(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
		}
		for (its=0;its<MAXIT;its++) { //refinement by Newton's method
			p1=1.0;
			p2=0.0;
			for (j=0;j<n;j++) { //Loop up the recurrence relation to get the Laguerre polynomial evaluated at z.
				p3=p2;
				p2=p1;
				p1=((2*j+1+alf-z)*p2-(j+alf)*p3)/(j+1);
			}
			//p1 is now the desired Laguerre polynomial. We next compute pp, its derivative,
			//by a standard relation involving also p2, the polynomial of one lower order.
			pp=(n*p1-(n+alf)*p2)/z;
			z1=z;
			z=z1-p1/pp; //Newton's formula
			if (fabs(z-z1) <= EPS) 
				break;
		}
		if (its >= MAXIT) 
			errorMsg::reportError("too many iterations in gaulag");
		x[i]=z;
		w[i] = -exp(gammln(alf+n)-gammln(MDOUBLE(n)))/(pp*n*p2);
	}
}


void GLaguer::GetPhylipLaguer(const int categs, MDOUBLE alpha, Vdouble & points, Vdouble & weights)
{
  /* calculate rates and probabilities to approximate Gamma distribution
     of rates with "categs" categories and shape parameter "alpha" using
     rates and weights from Generalized Laguerre quadrature */

	points.resize(categs, 0.0);
	weights.resize(categs, 0.0);
	long i;
	raterootarray lgroot; /* roots of GLaguerre polynomials */
	double f, x, xi, y;

	alpha = alpha - 1.0;
	lgroot[1][1] = 1.0+alpha;
	for (i = 2; i <= categs; i++)
	{
		cerr<<lgroot[i][1]<<"\t";
		lgr(i, alpha, lgroot);                   /* get roots for L^(a)_n */
		cerr<<lgroot[i][1]<<endl;
	}
	/* here get weights */
	/* Gamma weights are (1+a)(1+a/2) ... (1+a/n)*x_i/((n+1)^2 [L_{n+1}^a(x_i)]^2)  */
	f = 1;
	for (i = 1; i <= categs; i++)
		f *= (1.0+alpha/i);
	for (i = 1; i <= categs; i++) {
		xi = lgroot[categs][i];
		y = glaguerre(categs+1, alpha, xi);
		x = f*xi/((categs+1)*(categs+1)*y*y);
		points[i-1] = xi/(1.0+alpha);
		weights[i-1] = x;
	}
}


void GLaguer::lgr(long m, double alpha, raterootarray lgroot)
{ /* For use by initgammacat.  Get roots of m-th Generalized Laguerre
     polynomial, given roots of (m-1)-th, these are to be
     stored in lgroot[m][] */
	long i;
	double upper, lower, x, y;
	bool dwn;   /* is function declining in this interval? */

	if (m == 1) {
		lgroot[1][1] = 1.0+alpha;
	} else {
		dwn = true;
		for (i=1; i<=m; i++) {
			if (i < m) {
				if (i == 1)
					lower = 0.0;
				else
					lower = lgroot[m-1][i-1];
				upper = lgroot[m-1][i];
			} 
			else { /* i == m, must search above */
				lower = lgroot[m-1][i-1];
				x = lgroot[m-1][m-1];
				do {
					x = 2.0*x;
					y = glaguerre(m, alpha,x);
				}	while ((dwn && (y > 0.0)) || ((!dwn) && (y < 0.0)));
				upper = x;
			}
			while (upper-lower > 0.000000001) {
				x = (upper+lower)/2.0;
				if (glaguerre(m, alpha, x) > 0.0) {
					if (dwn)
						lower = x;
					else
						upper = x;
				} 
				else {
					if (dwn)
						upper = x;
					else
						lower = x;
				}			 
			}
			lgroot[m][i] = (lower+upper)/2.0;
			dwn = !dwn; // switch for next one 
		}
	}
} /* lgr */


double GLaguer::glaguerre(long m, double b, double x)
{ /* Generalized Laguerre polynomial computed recursively.
     For use by initgammacat */
	long i;
	double gln, glnm1, glnp1; /* L_n, L_(n-1), L_(n+1) */

	if (m == 0)
		return 1.0;
	else {
		if (m == 1)
		return 1.0 + b - x;
		else {
			gln = 1.0+b-x;
			glnm1 = 1.0;
			for (i=2; i <= m; i++) {
				glnp1 = ((2*(i-1)+b+1.0-x)*gln - (i-1+b)*glnm1)/i;
				glnm1 = gln;
				gln = glnp1;
			}
		return gln;
		}
	}
} /* glaguerre */
