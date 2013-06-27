// $Id: numRec.cpp 5990 2009-03-19 10:21:20Z privmane $

#include "numRec.h"
#include "matrixUtils.h"
#include <cassert>
#include <algorithm>

#ifndef VERBOS
#define VERBOS
#endif

void validateSym(VVdouble & v) {
	const MDOUBLE epsilon = 0.00000001;
	for (int i=0; i < v.size(); ++i) {
		for (int j=i+1; j < v.size(); ++j) {
			if (fabs(v[i][j] - v[j][i])> epsilon) {
				LOG(5,<<"v["<<i<<"]["<<j<<"]="<<v[i][j]<<endl);
				LOG(5,<<"v["<<j<<"]["<<i<<"]="<<v[j][i]<<endl);

				errorMsg::reportError("trying to find eigen values to non-sym matrix");
			}
			else v[i][j] = v[j][i];
		}
	}
}

int MyJacobi(VVdouble &Insym, VVdouble &RightEigenV, Vdouble &EigenValues) {
	validateSym(Insym);
	const int MaxNumberOfSweeps = 100000;
	VVdouble& v = RightEigenV;
	VVdouble& a = Insym;
	Vdouble& d = EigenValues;
	//CheckSizeAndTypeAndResizeIfNessary();
	int i,j;
	const int size = v.size();
		
	// preparing V to be the indentity matrix
	for (i=0; i<size; ++i) {
	  for (int j=0; j<size ; ++j) v[i][j]=0.0;
	  v[i][i] = 1.0;
	}
	
		
	for (i=0 ; i<size; ++i ) {
	  d[i] = a[i][i];
	}

	MDOUBLE sm = 0.0; // sm is the sum of the off-diagonal elements
	int ip, iq;
	for (i = 0; i< MaxNumberOfSweeps ; ++i) {
		sm = 0.0;
		for (ip = 0; ip<size ; ++ip) {
			for (iq = ip+1; iq <size; ++iq)	sm +=fabs (a[ip][iq]);
		}
		//if(i%300==0)
		//	LOG(5,<<"sm= "<<sm<<endl);
		if (sm == 0.0) return 0; // the program is suppose to return here, after some rounds of i.
		MDOUBLE tresh;
		if (i<3) tresh = 0.2 * sm / (size*size); else tresh = 0.0;

		MDOUBLE g;
		for (ip=0 ; ip<size; ++ip) {
			for (iq = ip+1 ; iq<size; ++iq) {
				g = 100.0*fabs(a[ip][iq]);

#ifdef VERBOS
	if (g<10e-50) {
		LOG(5,<<"small g!"<<endl);
		if ((i>3 && (fabs(d[ip]+g) == fabs(d[ip])) && (fabs(d[iq]+g)==fabs(d[iq])))==false) {
			LOG(5,<<"g is small: "<<g<< "yes, it is not zeroed"<<endl);
			LOG(5,<<"because d[ip] is: "<<d[ip]<<" and d[iq] is: "<<d[iq]<<endl);
			LOG(5,<<"ip is: "<<ip<<" iq is: "<<iq<<endl);
		}
	}
#endif //VERBOS
				if (i>3 && (fabs(d[ip]+g) == fabs(d[ip])) && (fabs(d[iq]+g)==fabs(d[iq])) ) {
					a[ip][iq] = 0.0;
				}
				else if (fabs(a[ip][iq]) > tresh) {
					MDOUBLE h;
					MDOUBLE t;
					MDOUBLE theta;
					h = d[iq]-d[ip];
	//  assert(h!=0);
					if (fabs(h) + g == fabs(h)) {
						assert(h!=0);
						t = a[ip][iq] / h;
					}
					else {
						theta = 0.5*h/(a[ip][iq]);
						t = 1.0 / (fabs(theta)+sqrt(1.0+theta*theta));
						if (theta<0.0) t = -t;
					}
					MDOUBLE c,s;
					c = 1.0 / sqrt(1.0+t*t);
					s = t*c;
					MDOUBLE tau;
					tau = s/ (1.0 + c);
					h = t * a[ip][iq];
					
					d[ip] = d[ip] - t * a[ip][iq];
					d[iq] = d[iq] + t * a[ip][iq];
					a[ip][iq]=0.0;
					MDOUBLE tmp1, tmp2;
					for (j = 0; j < ip; ++j) {
					  tmp1 = a[j][ip] - s*(a[j][iq]+a[j][ip]*tau); // updating the above element of a...
					  tmp2 = a[j][iq] + s*(a[j][ip]-a[j][iq]*tau);
					  a[j][ip] = tmp1; 
					  a[j][iq] = tmp2;
					}
									
					for (j = ip+1;j<iq; ++j) {
					  tmp1 = a[ip][j] - s*(a[j][iq]+a[ip][j]*tau); // updating the above element of a..
					  tmp2 = a[j][iq] + s*(a[ip][j]-a[j][iq]*tau);
					  a[ip][j] = tmp1;
					  a[j][iq] = tmp2;
					}
					
					for (j = iq+1; j< size ; ++j) {
					  tmp1 = a[ip][j] - s*(a[iq][j]+a[ip][j]*tau); // updating the above element of a..
					  tmp2 = a[iq][j] + s*(a[ip][j]-a[iq][j]*tau);
					  a[ip][j] = tmp1;
					  a[iq][j] = tmp2;
					}
									
					for (j = 0; j< size ; ++j) {
					  tmp1 = v[j][ip] - s*(v[j][iq]+v[j][ip]*tau); // updating v
					  tmp2 = v[j][iq] + s*(v[j][ip]-v[j][iq]*tau);
					  v[j][ip] = tmp1;
					  v[j][iq] = tmp2;
					}
				} // end of "else if (fabs(a[ip][iq] > tresh)"
			} // end of for (iq = ...
		} // end of for (ip = ...
	} // end of for (i = 0; i< MaxNumberOfSweeps ; ++i) {
	vector<string> err;
	err.push_back("problems in function MyJacobi. more than MaxNumberOfSweeps were necesary.");
	errorMsg::reportError(err);
	
	return -1;
} //end of function





///////////////////////////////////////////
//Adi cahnges	//////////////////////////
/////////////////////////////////////////
MDOUBLE sign(MDOUBLE a,MDOUBLE b){ 
	return (b>0?fabs(a):-fabs(a));
}

MDOUBLE pythag(const MDOUBLE a, const MDOUBLE b){
	return sqrt(pow(a,2)+pow(b,2));
}


void houseHolder(VVdouble &mat,VVdouble &Q){
	MDOUBLE sigma=0,H,sqrtSigma,K=0,tmp;
	int c,r,j,i,n = mat.size();
	Q.resize(n);
	for(i=0;i<n;i++){
		Q.resize(n);
	}
	for (i=0;i<n;i++)
		Q[i].resize(n,0.0);
	Vdouble p,q,u;
	p.resize(n,0.0);
	q.resize(n,0.0);
	u.resize(n,0.0);
	for (i=n-1;i>1;i--){
		sigma=0; //init sigma
		K=0;	//init K
	
		for(j=0;j<i;j++)
			sigma+= mat[i][j]*mat[i][j]; //compute sigma: O(n)

		sqrtSigma = mat[i][i-1]>=0.0 ? sqrt(sigma) : -sqrt(sigma); //compute sqrt of sigma +/-

		H=sigma+mat[i][i-1]*sqrtSigma; //comute H = 0.5*|u|^2.  until here O(n)
	
		/***createing U*******/
		for(r=0;r<i;r++) {   //update vector u with row i the matrix until i; //takes  O(n)
			Q[i][r]= u[r] = mat[i][r];
			Q[r][i] = u[r]/H;
		}
		u[i-1]+=sqrtSigma; //update element (i,i-1)
		Q[i][i-1]=u[i-1];
		Q[i-1][i]=u[i-1]/H;
		for(r=i;r<n;r++) //update elemnts (i,j) =0 for j>=i.
			u[r]=0.0;
		/***********************/
		for(r=0;r<n;r++){ //compute vector p O(n^2)
			p[r]=0.0;
			for (c=0;c<i;c++)
				p[r]+=mat[r][c]*u[c]; //compute AU		
			p[r]/=H; // ->AU/H		
		}
	
		for(r=0;r<i;r++) // compure K O(n)
			K+=u[r]*p[r];
		K/=(2*H);
	//	cout<<"K is: "<<K<<endl;
		
		for(r=0;r<n;r++) //compute vector q O(n)
			q[r]=p[r]-K*u[r];
		
		for(r=0;r<=i;r++) {//update matrix O(n^2) only part of the matrix
			for(c=0;c<=i;c++)
				mat[r][c]-=q[r]*u[c]+u[r]*q[c];
		}
					
	}
	for (i=0;i<n;i++){
		for(j=0;j<i;j++){
			tmp=0;
			for(c=0;c<i;c++)
				tmp+=Q[i][c]*Q[c][j];
			for(c=0;c<i;c++)
				Q[c][j]-=tmp*Q[c][i];
		}
		Q[i][i]=1;
		for(j=0;j<i;j++)
			Q[j][i]=Q[i][j]=0.0;
	}
}

void tred2(VVdouble &a, Vdouble &d, Vdouble &e) //a = symmetricMatrix,d = diagonal,e = offdiagonal
{
	int l,k,j,i;
	MDOUBLE scale,hh,h,g,f;

	int n=d.size();
	for (i=n-1;i>0;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<l+1;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=0;k<l+1;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=0;j<l+1;j++) {
				// Next statement can be omitted if eigenvectors not wanted
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=0;k<j+1;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<l+1;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=0;j<l+1;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<j+1;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	// Next statement can be omitted if eigenvectors not wanted
	d[0]=0.0;
	e[0]=0.0;
	// Contents of this loop can be omitted if eigenvectors not
	//	wanted except for statement d[i]=a[i][i];
	for (i=0;i<n;i++) {
		l=i;
		if (d[i] != 0.0) {
			for (j=0;j<l;j++) {
				g=0.0;
				for (k=0;k<l;k++)
					g += a[i][k]*a[k][j];
				for (k=0;k<l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=0;j<l;j++) a[j][i]=a[i][j]=0.0;
	}
}

//called if houseHolder was used - the modified QL implementation corresponding to the modified implementation of householder
/*
void QL(Vdouble &d, Vdouble &e, VVdouble &z){
	int m,l,iter,i,k;
	MDOUBLE s,r,p,g,f,dd,c,b;
		
	int n=d.size();
//*	for (i=1;i<n;i++) e[i-1]=e[i];
//*	e[n-1]=0.0;
//*	e.push_back(0);//since in my algorithm I return an n-1 sized e
	for (l=0;l<n;l++) {
		iter=0;
		do {
			for (m=l;m<n-1;m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if (fabs(e[m])+dd == dd) break;
			}
			if (m != l) {
				if (iter++ == 30) errorMsg::reportError("Too many iterations in QL");
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+sign(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					// Next loop can be omitted if eigenvectors not wanted
					for (k=0;k<n;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}
*/


//called if tred2 was used - the original QL implementation from numerical recepies
void QL(Vdouble &d, Vdouble &e, VVdouble &z){
	int m,l,iter,i,k;
	MDOUBLE s,r,p,g,f,dd,c,b;

	int n=d.size();
	for(i=1;i<n;i++){
		e[i-1]=e[i];
	}
	e[n-1]=0.0;
	for(l=0;l<n;l++){
		iter=0;
		do {
			for(m=l;m<n-1;m++){
				dd=fabs(d[m])+fabs(d[m+1]);
				if(fabs(e[m])+dd == dd) break;
			}
			if(m!=l){
				if(iter++==30){
					errorMsg::reportError("too many iteration in QL");
				}
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+sign(r,g));
				s=c=1.0;
				p=0.0;
				for(i=m-1;i>=l;i--){
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if(r==0.0){
						d[i+1]-=p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					for(k=0;k<n;k++){
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if(r==0 && i>=l) continue;
				d[l]-=p;
				e[l]=g;
				e[m]=0.0;
			}
		}
		while(m!=l);
	}
}



/************************************************************************/
//diaganol will be eigen values and fill matrix of eigen vectors.                                                                                   */
/************************************************************************/

//A modified implementation for eigen analysis, using the house holder function.
/*
void computeEigenSystem(VVdouble &symmetricMatrix,VVdouble &eigenVectros,Vdouble &diagonal){

	houseHolder(symmetricMatrix,eigenVectros);

	Vdouble offdiagonal;
	offdiagonal.resize(symmetricMatrix.size());
	for (int i=0; i<symmetricMatrix.size(); i++){
		diagonal[i]=symmetricMatrix[i][i];
	}
	for (int i2=0; i2<symmetricMatrix.size()-1; i2++){
		offdiagonal[i2]=symmetricMatrix[i2+1][i2];
	}
	
	QL(diagonal,offdiagonal,eigenVectros);
	return;
}
*/

//Uses original implementation of tred2 function for eigen analysis, copied from numerical recepies p474. 
void computeEigenSystem(VVdouble &symmetricMatrix,VVdouble &eigenVectros,Vdouble &diagonal){
	
	Vdouble offdiagonal;
	offdiagonal.resize(symmetricMatrix.size());

	tred2(symmetricMatrix,diagonal,offdiagonal);

	eigenVectros = symmetricMatrix;

	QL(diagonal,offdiagonal,eigenVectros);

	return;
}


// the following two functions used for Kolomogorov-Smirnoff test
MDOUBLE performKSTest(const uniformDistribution& empiricalDist, Vdouble& observedDist)
{
	MDOUBLE pVal = 0.0;
	MDOUBLE distance = 0.0;

	int j;
	MDOUBLE dt,en,fn,fo = 0.0;

	int n = observedDist.size();
	sort(observedDist.begin(),observedDist.end());
	en = n;
	MDOUBLE cdfObserved = 0.0;
	for(j = 0; j < n; ++j){
		cdfObserved+=observedDist[j];
		fn = (j+1)/en;
		dt = max(fabs(fo-cdfObserved),fabs(fn-cdfObserved));
		if(dt > distance)
			distance = dt;
		fo = fn;
	}
	en = sqrt(en);
	pVal = computeProbForKS((en+0.12+0.11/en)*distance);
	return pVal;
}

// function called only by performKSTest
MDOUBLE computeProbForKS (const MDOUBLE QsParam)
{
	const MDOUBLE EPS1 = 1.0e-6,EPS2 = 1.0e-16;
	int j;
	MDOUBLE a2,fac = 2.0, sum = 0.0, term, termbf = 0.0;

	a2 = -2.0*QsParam*QsParam;
	for(j = 1; j <= 100; ++j){
		term = fac*exp(a2*j*j);
		sum += term;
		if(fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum)
			return sum;
		fac = -fac;
		termbf = fabs(term);
	}
	return 1.0; //get here only by failing to converge
}
