// $Id: readDatMatrix.cpp 5805 2009-01-20 09:19:26Z adido $
//#ifndef unix
//#define SSTREAM_KNOWN
//#endif

//#ifdef SSTREAM_KNOWN
#include <sstream>
//#else 
//#include <strstream> //oldVersion
//#endif


#include <cassert>
#include "readDatMatrix.h"
#include "errorMsg.h"
#include "logFile.h"

//#define VERBOS

void normalizeQ(VVdouble& q, const Vdouble& freq) {
	MDOUBLE sum =0;
	int i=0,j=0;
	for (i=0; i < q.size(); ++i) {
		sum += q[i][i]*freq[i];
	}
	assert(sum!=0);
	MDOUBLE oneDividedBySum = -1.0/sum; // to avoid many divisions.

	for (i=0; i < q.size(); ++i) {
		for (j=0; j < q.size(); ++j) {
			q[i][j] = q[i][j]*oneDividedBySum;
		}
	}
}

void readDatMatrixFromFile(const string & matrixFileName,
						   VVdouble & subMatrix,
						   Vdouble & freq) {
	cout<<"****readDatMatrixFromFile******"<<endl;
	int i=0,j=0; //indices
	ifstream in(matrixFileName.c_str());
	if (!in) {
		errorMsg::reportError("unable to open matrix data file");
	}
	
	int alphaSize;
	if (matrixFileName == "adrianCodon.dat.q")
		alphaSize = 61;
	else 
		alphaSize = 20;
	subMatrix.resize(alphaSize);
	for ( i=0; i < alphaSize; ++i) subMatrix[i].resize(alphaSize,0.0);
	freq.resize(alphaSize,0.0);

	for (i=1; i < subMatrix.size(); ++i) {
		for (j=0; j <i;++j) {
			in>>subMatrix[i][j];
			subMatrix[j][i] = subMatrix[i][j];
		}
	}
	for (i=0; i < subMatrix.size(); ++i) {
		in>>freq[i];
	}
	in.close();

	//check:
	//LOG(5,<<" priting the 5*5 top part of the sub matrix: "<<endl);
	//for (i=0; i < 5; ++i) {
	//	for (j=0; j <5;++j) {
	//		LOG(5,<<subMatrix[i][j]<<" ");
	//	}
	//	LOG(5,<<endl);
	//}
	//LOG(5,<<"the 5 last freqs: "<<endl);
	//for (i=15; i < 20; ++i) {
	//	LOG(5,<<freq[i]<<" ");
	//}
}

void readDatMatrixFromString(const string & matrixFileString,
			     VVdouble & subMatrix,
			     Vdouble & freq, int alphaSize) {
	int i=0,j=0; //indices	
	//#ifdef SSTREAM_KNOWN
	stringstream in(matrixFileString.c_str());
// #else
// 	istrstream in(matrixFileString.c_str()); // OLD VERSION
//#endif
	if (!in) {
		errorMsg::reportError("unable to open matrix data buffer");
	}

	
	subMatrix.resize(alphaSize);
	for ( i=0; i < alphaSize; ++i) subMatrix[i].resize(alphaSize,0.0);
	freq.resize(alphaSize,0.0);

	for (i=1; i < alphaSize; ++i) {
		for (j=0; j <i;++j) {
			in>>subMatrix[i][j];
			subMatrix[j][i] = subMatrix[i][j];
		}
	}
	for (i=0; i < alphaSize; ++i) {
		in>>freq[i];
	}
}


#include "fromQtoPt.h"
#include "definitions.h"

#include <iostream>
using namespace std;

void pupAll::fillMatricesFromFile(const string & dataFileString) {
	VVdouble sMatrix;
	readDatMatrixFromFile(dataFileString,sMatrix,_freq);
	//	readDatMatrixFromString(dataFileString,sMatrix,_freq);
	VVdouble qMatrix = fromWagSandFreqToQ(sMatrix,_freq);
	
	q2pt q2pt1;
	q2pt1.fillFromRateMatrix(_freq,qMatrix);
	_leftEigen = q2pt1.getLeftEigen();
	_rightEigen = q2pt1.getRightEigen();
	_eigenVector = q2pt1.getEigenVec();
}
void pupAll::fillMatricesFromFile(const string & dataFileString, const Vdouble & freq) {
#ifdef VERBOS
	LOG(5,<<"dataFileString = "<<dataFileString<<endl);
#endif

	VVdouble sMatrix;
	readDatMatrixFromFile(dataFileString,sMatrix,_freq);
	_freq=freq;
	VVdouble qMatrix = fromWagSandFreqToQ(sMatrix,_freq);
	
	q2pt q2pt1;
	q2pt1.fillFromRateMatrix(_freq,qMatrix);
	_leftEigen = q2pt1.getLeftEigen();
	_rightEigen = q2pt1.getRightEigen();
	_eigenVector = q2pt1.getEigenVec();
}

void pupAll::fillMatrices(const string & dataFileString,int alphaSize) {
	VVdouble sMatrix;
	readDatMatrixFromString(dataFileString,sMatrix,_freq,alphaSize);
	//	readDatMatrixFromString(dataFileString,sMatrix,_freq);
	VVdouble qMatrix = fromWagSandFreqToQ(sMatrix,_freq);
	
	q2pt q2pt1;
	q2pt1.fillFromRateMatrix(_freq,qMatrix);
	_leftEigen = q2pt1.getLeftEigen();
	_rightEigen = q2pt1.getRightEigen();
	_eigenVector = q2pt1.getEigenVec();
}
void pupAll::fillMatrices(const string & dataFileString, const Vdouble & freq) {
	VVdouble sMatrix;
	readDatMatrixFromString(dataFileString,sMatrix,_freq);
	_freq=freq;
	VVdouble qMatrix = fromWagSandFreqToQ(sMatrix,_freq);
	
	q2pt q2pt1;
	q2pt1.fillFromRateMatrix(_freq,qMatrix);
	_leftEigen = q2pt1.getLeftEigen();
	_rightEigen = q2pt1.getRightEigen();
	_eigenVector = q2pt1.getEigenVec();
}

const MDOUBLE pupAll::Pij_t(const int i, const int j, const MDOUBLE t) const {
	if (t<0) {
		LOG(5,<<"negative length in routine Pij_t "<<endl);
		LOG(5,<<" t = " <<t<<endl);
		errorMsg::reportError("negative length in routine Pij_t");
	}
//	if ((_freq[i] == 0.0) || (_freq[j] == 0.0)) return 0.0;
	MDOUBLE sum=0;
	int alphaSize = _freq.size();
	for (int k=0 ; k<alphaSize ; ++k) {
		sum+=( _leftEigen[i][k]*_rightEigen[k][j]*exp(_eigenVector[k]*t) );
	}
	if (currectFloatingPointProblems(sum)) return sum; 
//	LOG(1,<<"err Pij_t i="<<i<<" j= "<<j<<" dis= "<<t<<" res= "<<sum<<endl);//sum is not in [0,1]
	errorMsg::reportError("error in function pijt... ");return 0;
}

const MDOUBLE pupAll::dPij_dt(const int i,const  int j, const MDOUBLE t) const {
//	if ((_freq[i] == 0.0) || (_freq[j] == 0.0)) return 0.0;
	MDOUBLE sum=0;
	int alphaSize = _freq.size();
	for (int k=0 ; k<alphaSize ; ++k) {
		sum+=( _leftEigen[i][k]*_rightEigen[k][j]*exp(_eigenVector[k]*t)*_eigenVector[k]);
	}
	return sum;
}


const MDOUBLE pupAll::d2Pij_dt2(const int i,const int j, const MDOUBLE t) const {
//	if ((_freq[i] == 0.0) || (_freq[j] == 0.0)) return 0.0;
	MDOUBLE sum=0;;
	int alphaSize = _freq.size();
	for (int k=0 ; k<alphaSize ; ++k) {
		sum+=( _leftEigen[i][k]*_rightEigen[k][j]*exp(_eigenVector[k]*t)*_eigenVector[k]*_eigenVector[k]);
	}
	return sum;
}
// this gives the likelihood of j given i at distance t and gamma
// parameter alpha.  The result presented here is the integral over the
// rates (according to the gamma distribution with parameter alpah).  see Yang's (93) paper.
const MDOUBLE pupAll::Pij_tAlpha(const int i, const int j, const MDOUBLE t, const MDOUBLE alpha) const {
	if (t<0) {
		LOG(5,<<"negative length in routine Pij_tAlpha "<<endl);
		LOG(5,<<" t = " <<t<<endl);
		errorMsg::reportError("negative length in routine Pij_tAlpha");
	}
	MDOUBLE sum=0;
	for (int k=0 ; k<20 ; ++k) {
		sum+=( _leftEigen[i][k]*_rightEigen[k][j]*pow(1-_eigenVector[k]*t/alpha,-alpha));
	}
	if (currectFloatingPointProblems(sum)) return sum; 
	errorMsg::reportError("error in function pijtAlpha... ");return 0;
}


const MDOUBLE pupAll::Pij_tAlpha_dt(const int i, const int j, const MDOUBLE t, const MDOUBLE alpha) const {
	if (t<0) {
		LOG(5,<<"negative length in routine Pij_tAlpha_dt "<<endl);
		LOG(5,<<" t = " <<t<<endl);
		errorMsg::reportError("negative length in routine Pij_tAlpha_dt");
	}
	MDOUBLE sum=0;
	for (int k=0 ; k<20 ; ++k) {
		sum+=( _leftEigen[i][k]*_rightEigen[k][j]* _eigenVector[k]*  pow(1-_eigenVector[k]*t/alpha,-alpha-1));
	}
	return sum; 
}
const MDOUBLE pupAll::Pij_tAlpha_dt2(const int i, const int j, const MDOUBLE t, const MDOUBLE alpha) const {
	if (t<0) {
		LOG(5,<<"negative length in routine Pij_tAlpha_dt2 "<<endl);
		LOG(5,<<" t = " <<t<<endl);
		errorMsg::reportError("negative length in routine Pij_tAlpha_dt2");
	}
	MDOUBLE sum=0;
	for (int k=0 ; k<20 ; ++k) {
		sum+=( _leftEigen[i][k]*_rightEigen[k][j]* (1+1/alpha) *_eigenVector[k]*_eigenVector[k]*  pow(1-_eigenVector[k]*t/alpha,-alpha-2));
	}
	return sum; 
}

bool pupAll::currectFloatingPointProblems(MDOUBLE& sum) const {
	if ((sum * (sum+err_allow_for_pijt_function))<0) sum=0;
	if (((sum-1) * (sum-1.0-err_allow_for_pijt_function))<0) sum=1;
	if ((sum>1) || (sum<0)) return false;
	return true;
}

VVdouble fromWagSandFreqToQ(const VVdouble & s,const Vdouble& freq){
	VVdouble q(s.size());
	for (int z=0; z < q.size(); ++z) q[z].resize(s.size(),0.0);
	int i,j;
	MDOUBLE sum;
	for ( i=0; i < s.size(); ++i) {
		sum =0;
		for (j=0; j < s.size(); ++j) {
			if (i!=j) q[i][j] = s[i][j]* freq[j];
			sum += q[i][j];
		}
		q[i][i] = -sum;
	}

	// normalizing q:
	normalizeQ(q,freq);


	// check:
	//sum =0;
	//for (i=0; i < s.size(); ++i){
	//	sum += q[i][i]*freq[i];
	//}
	//LOG(5,<<" SUM OF DIAGOPNAL Q IS (should be -1) "<<sum<<endl);
	return q;

}

