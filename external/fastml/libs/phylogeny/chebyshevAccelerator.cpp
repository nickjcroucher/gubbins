// $Id: chebyshevAccelerator.cpp 962 2006-11-07 15:13:34Z privmane $

#include "chebyshevAccelerator.h"
#include <cmath>
#include <cassert>

chebyshevAccelerator::chebyshevAccelerator(const chebyshevAccelerator& other):
	_alphabetSize(other._alphabetSize),
	_totalNumOfCoef(other._totalNumOfCoef),
	_usingNumberOfCoef(other._usingNumberOfCoef),
	_pb(NULL),
	_rightRange(other._rightRange),
	_leftRange(other._leftRange){
		if (other._pb != NULL) _pb = other._pb->clone();
		chebi_coff=other.chebi_coff;
		chebi_dervation_coff=other.chebi_dervation_coff;
		chebi_sec_dervation_coff=other.chebi_sec_dervation_coff;
}

chebyshevAccelerator::chebyshevAccelerator(
	 replacementModel* pb,
	const int alphanetSize,
	const int totalNumOfCoef,
	const int usingNumberOfCoef,
	const MDOUBLE rightRange,
	const MDOUBLE leftRange
	): _alphabetSize(alphanetSize),
	_totalNumOfCoef(totalNumOfCoef), _usingNumberOfCoef(usingNumberOfCoef),_pb(pb->clone()), _rightRange(rightRange), _leftRange(leftRange)
//----------------------------------------------------------------------------------
//input:	non
//output:	non
//doing:	filling the member chebi_coff[][][]; chebi_coff[1][2][4] is the forth 
//			chebichev coefficient in the chebichev polynom of the function 
//			slow_pij(1,2,t);			
//----------------------------------------------------------------------------------
{
	int tmp, tmp1;
	for (tmp  = 0; tmp < _alphabetSize ; tmp ++) {
		
		chebi_coff.resize(_alphabetSize);
		chebi_dervation_coff.resize(_alphabetSize);
		chebi_sec_dervation_coff.resize(_alphabetSize);

		for (tmp1  = 0; tmp1 < _alphabetSize ; tmp1 ++) {
			chebi_coff[tmp].resize(_alphabetSize);
			chebi_dervation_coff[tmp].resize(_alphabetSize);
			chebi_sec_dervation_coff[tmp].resize(_alphabetSize); 
			for (tmp1  = 0; tmp1 < _alphabetSize ; tmp1 ++) {
				chebi_coff[tmp][tmp1].resize(_totalNumOfCoef);
				chebi_dervation_coff[tmp][tmp1].resize(_totalNumOfCoef);
				chebi_sec_dervation_coff[tmp][tmp1].resize(_totalNumOfCoef); 
			}
		}
	}


	Vdouble coffij(_totalNumOfCoef);
	Vdouble coffij_of_derviation(_totalNumOfCoef);
	Vdouble coffij_of_second_derivation(_totalNumOfCoef);

	
	for (int from_aa =0; from_aa<_alphabetSize ; ++ from_aa)
	{
		for (int to_aa =0; to_aa<_alphabetSize ; ++ to_aa)
		{
			chebft(coffij,_totalNumOfCoef,from_aa,to_aa);
			chder(coffij,coffij_of_derviation,_totalNumOfCoef);
			chder(coffij_of_derviation,coffij_of_second_derivation,_totalNumOfCoef);

			for (int tmp=0; tmp<_totalNumOfCoef;++tmp)
			{
				chebi_coff[from_aa][to_aa][tmp] = coffij[tmp];
				chebi_dervation_coff[from_aa][to_aa][tmp] = coffij_of_derviation[tmp];
				chebi_sec_dervation_coff[from_aa][to_aa][tmp] = coffij_of_second_derivation[tmp];
			}

		}
	}
}


void chebyshevAccelerator::chebft(Vdouble& c, int n, int from_aa, int to_aa) {
//----------------------------------------------------------------------------------
//input:	c[] is the vector where the cofficient will be
//			from aa and to_aa are for chosing the right function to be developed
//output:	non
//doing:	calculating the  chebichev coefficient in the chebichev polynom of the function 
//			slow_pij(from_aa,to_aa,t), and put them in the c[] vector			
//----------------------------------------------------------------------------------
	int k,j;
	MDOUBLE fac,bpa,bma;

	Vdouble f;
	f.resize(n);
	bma=0.5*(_rightRange-_leftRange);
	bpa=0.5*(_rightRange+_leftRange);
	for (k=0;k<n;k++) {
		MDOUBLE y=cos(3.141592653589793*(k+0.5)/n);
		f[k]= _pb->Pij_t(from_aa,to_aa,y*bma+bpa); //(*func)(y*bma+bpa);
	}
	fac=2.0/n;
	for (j=0;j<n;j++) {
		MDOUBLE sum=0.0;
		for (k=0;k<n;k++)
			sum += f[k]*cos(3.141592653589793*j*(k+0.5)/n);
		c[j]=fac*sum;
	}
	
}


const MDOUBLE chebyshevAccelerator::Pij_t(const int from_aa, const int to_aa, const MDOUBLE x) const
//----------------------------------------------------------------------------------
//input:	like pij_t
//output:	the probabilty
//doing:	calculating with the polinom of chebi and via eigenvalue decomposition			
//----------------------------------------------------------------------------------
{
  
	MDOUBLE d=0.0,dd=0.0,sv,y,y2,check;
	int j;

	if ((x-_leftRange)*(x-_rightRange) > 0.0) {
	return _pb->Pij_t(from_aa,to_aa,x);
//		errorMsg::reportError("x not in range in routine fast_Pij_t");// also quit the program
	}

	y2=2.0*(y=(2.0*x-_leftRange-_rightRange)/(_rightRange-_leftRange));
	for (j=_usingNumberOfCoef;j>0;j--) {
		sv=d;
		d=y2*d-dd+chebi_coff[from_aa][to_aa][j];
		dd=sv;
	}
	check =  y*d-dd+0.5*chebi_coff[from_aa][to_aa][0];
	if ((check>1) || (check<=0)) check = _pb->Pij_t(from_aa,to_aa,x);
	assert(check<=1);
	assert(check>=0);
	return check;
}


const MDOUBLE chebyshevAccelerator::dPij_dt(const int from_aa, const int to_aa, const MDOUBLE x) const
//----------------------------------------------------------------------------------
//input:	like pij_t
//output:	the derivation of probabilty
//doing:	calculating with the polinom of chebi and via eigenvalue decomposition			
//----------------------------------------------------------------------------------
{
  
	MDOUBLE d=0.0,dd=0.0,sv,y,y2;
	int j;

	if ((x-_leftRange)*(x-_rightRange) > 0.0) {
		return _pb->dPij_dt(from_aa,to_aa,x);
	}
	y2=2.0*(y=(2.0*x-_leftRange-_rightRange)/(_rightRange-_leftRange));
	for (j=_usingNumberOfCoef;j>0;j--) {
		sv=d;
		d=y2*d-dd+chebi_dervation_coff[from_aa][to_aa][j];
		dd=sv;
	}
	return y*d-dd+0.5*chebi_dervation_coff[from_aa][to_aa][0];
}


const MDOUBLE chebyshevAccelerator::d2Pij_dt2(const int from_aa, const int to_aa, const MDOUBLE x) const {
//----------------------------------------------------------------------------------
//input:	like pij_t
//output:	the second derivation of the probabilty
//doing:	calculating with the polynom of chebi and via eigenvalue decomposition			
//----------------------------------------------------------------------------------
	MDOUBLE d=0.0,dd=0.0,sv,y,y2;
	int j;

	if ((x-_leftRange)*(x-_rightRange) > 0.0) {
			return _pb->d2Pij_dt2(from_aa,to_aa,x);
	}
	y2=2.0*(y=(2.0*x-_leftRange-_rightRange)/(_rightRange-_leftRange));
	for (j=_usingNumberOfCoef;j>0;j--) {
		sv=d;
		d=y2*d-dd+chebi_sec_dervation_coff[from_aa][to_aa][j];
		dd=sv;
	}
	return y*d-dd+0.5*chebi_sec_dervation_coff[from_aa][to_aa][0];
}




void chebyshevAccelerator::chder(Vdouble &c, Vdouble &cder, int n) {
//----------------------------------------------------------------------------------
//input:	chebicev coff of f(x) i.e. in c[]. n is the vector size
//output:	chebicev coff of df(x)/dx i.e. in cder[]
//doing:	calculating the coff of the dervation from the coff of f.
//reference:numercal recepies in c, pg 195.			
//----------------------------------------------------------------------------------
	int j;
	MDOUBLE con;

	cder[n-1]=0.0;
	cder[n-2]=2*(n-1)*c[n-1];
	for (j=n-3;j>=0;j--)
		cder[j]=cder[j+2]+2*(j+1)*c[j+1];
	con=2.0f/(_rightRange-_leftRange);
	for (j=0;j<n;j++)
		cder[j] *= con;
}





