// $Id: fromQtoPt.cpp 5788 2009-01-19 22:24:16Z rubi $

#include "definitions.h"
#include "fromQtoPt.h"
#include "errorMsg.h"
#include "numRec.h"
#include "matrixUtils.h"
#include <iostream>
using namespace std;
#include <cassert>

//#define VERBOS




void q2pt::fillFromRateMatrix(const vector<MDOUBLE>& freq,
		   const VVdouble & qMatrix) {
	   // we first decompose Q to (F^0.5) M (F^-0.5)
	   // F is a diagonal matrix of the frequencies
	   // M is the symetrical matrix representation of Q.
	
	VVdouble q_sym;
	const int matrix_size = qMatrix.size();
	q_sym.resize(matrix_size);
	int k=0;
	for (k=0; k < q_sym.size(); ++k) q_sym[k].resize(matrix_size);
	calc_symmetric_q(qMatrix,q_sym,freq);
	// now we have to find the eigen-vector decomposition of the q_sym.
	VVdouble v; // v is the eigen vectors of the symetrical matrix.
	v.resize(matrix_size);
	for (k=0; k < qMatrix.size(); ++k) v[k].resize(matrix_size);
	Vdouble eigenValues(matrix_size);
	
	// symmetric_1pam = [v] [eigenValues] [transpose(v)]
	//MyJacobi(q_sym,v, eigenValues); // notice that inv([v]) = [v] transpose;
	

	/////i changed
	computeEigenSystem(q_sym,v,eigenValues);

	////
//#ifdef VERBOS
//	LOG(5,<<"The eigen-vector matrix of the decomposition of the symetric matrix\n");
//	for (int k1=0; k1 < v.size(); ++k1) {
//		for (int k2=0; k2<v[k1].size(); ++k2) {
//			LOG(5,<<v[k1][k2]<<" ");
//		}
//		LOG(5,<<endl);
//	}
//#endif 


	VVdouble left_eig_of_pam; // v is the eigen vectors of the symetrical matrix.
	left_eig_of_pam.resize(matrix_size);
	for (k=0; k < left_eig_of_pam.size(); ++k) left_eig_of_pam[k].resize(matrix_size);
	VVdouble right_eig_of_pam; // v is the eigen vectors of the symetrical matrix.
	right_eig_of_pam.resize(matrix_size);
	for (k=0; k < right_eig_of_pam.size(); ++k) right_eig_of_pam[k].resize(matrix_size);

	calc_left_and_right_eig_of_pam(left_eig_of_pam,right_eig_of_pam,v,freq);
	
	_leftEigen=left_eig_of_pam;
	_rightEigen=right_eig_of_pam;
	_eigenVector=eigenValues;
	Vdouble _freq=freq;
	// printing a pij(1);
	//MDOUBLE t = 1;
	//string fileName = "D://My Documents//adid//nimrod//inputs//inputs//aligned tce//aligned tce//P.F//P.F. vs P.F//eigenValues1.txt";
//	ofstream out(fileName.c_str());
//	for (int i=0;i<eigenValues.size();i++)
//		out<<eigenValues[i] <<" ";
//	out<<endl;
	//for (int aa1=0; aa1 < eigenValues.size(); ++aa1) {
	//	for (int aa2=0; aa2 < eigenValues.size(); ++aa2) {
	///		MDOUBLE sum=0;
	//		for (int k=0 ; k<eigenValues.size() ; ++k) {
	//			sum+=( left_eig_of_pam[aa1][k]*right_eig_of_pam[k][aa2]*exp(eigenValues[k]*t) );
	//		}
	//		LOG(5,<<sum<<" ");
//		}
//		LOG(5,<<endl);
//	}
}

void q2pt::fillFrom1PAMMatrix(const vector<MDOUBLE>& freq,const VVdouble & onePam)
{
	fillFromRateMatrix(freq,onePam);
	for (int i=0; i < 	_eigenVector.size(); ++i) {
		assert(_eigenVector[i]>0);
		_eigenVector[i] = log(_eigenVector[i])* 100;
	}
}

bool q2pt::currectFloatingPointProblems(MDOUBLE& sum) const {
	if ((sum * (sum+err_allow_for_pijt_function))<0) sum=0;
	if (((sum-1) * (sum-1.0-err_allow_for_pijt_function))<0) sum=1;
	if (!((sum<=1) && (sum>=0))) 
		return false;
	return true;
}

// Pij(t) = Sigma[k]{ [V]ik * [V^-1]kj * e^(Lamda_k*t) }
const MDOUBLE q2pt::Pij_t(const int i, const int j, const MDOUBLE t) const {
	if (t<0) errorMsg::reportError("negative length in routine Pij_t");
//	if ((_freq[i] == 0.0) || (_freq[j] == 0.0)) return 0.0;
	MDOUBLE sum=0;
	for (int k=0 ; k<_eigenVector.size() ; ++k) {
		sum+=( _leftEigen[i][k]*_rightEigen[k][j]*exp(_eigenVector[k]*t) );
	}
	if (currectFloatingPointProblems(sum)) return sum; 
//	LOG(1,<<"err Pij_t i="<<i<<" j= "<<j<<" dis= "<<t<<" res= "<<sum<<endl);//sum is not in [0,1]
	errorMsg::reportError("q2pt::Pij_t error in function pijt... ");return 0;
}

const MDOUBLE q2pt::dPij_dt(const int i,const  int j, const MDOUBLE t) const {
	MDOUBLE sum=0;
	for (int k=0 ; k<_eigenVector.size() ; ++k) {
		sum+=( _leftEigen[i][k]*_rightEigen[k][j]*exp(_eigenVector[k]*t)*_eigenVector[k]);
	}
	return sum;
}


const MDOUBLE q2pt::d2Pij_dt2(const int i,const int j, const MDOUBLE t) const {
	MDOUBLE sum=0;;
	for (int k=0 ; k<_eigenVector.size() ; ++k) {
		sum+=( _leftEigen[i][k]*_rightEigen[k][j]*exp(_eigenVector[k]*t)*_eigenVector[k]*_eigenVector[k]);
	}
	return sum;
}

void q2pt::calc_symmetric_q(const VVdouble &q_matrix,
							VVdouble &symmetric_q,
							const Vdouble & freq)
//----------------------------------------------------------------------------------
//input:	symmetric_1pam matrix is the output, pam1 is the input
//output:	non
//doing:	procedures to find eigen values work on symetrical matrices.
//			dayhoff 1 pam in a new basis is symetrical
//			the transformation is
//			
//			(1)  [symmetric_1pam] = [sqrt(pi)] * [pam1] * [1/sqrt(pi)]
//
//			[] for matrix. [sqrt(pi)] is a diagonal matrix were a[i][i] is the root of freq[i]
//reference: JME (1997) 45:696-703 Estimation of reversible substitution matrices from
//			 multiple pairs of sequences. Lars Arvestad and William J. Bruno.
//----------------------------------------------------------------------------------
{	
	int i,j;
	for (i=0; i<q_matrix.size(); ++i) {
		for (j=0; j<q_matrix.size(); ++j) {
			if (q_matrix[i][j] != 0.0) {
				 symmetric_q[i][j] = q_matrix[i][j]*sqrt(freq[i])/sqrt(freq[j]);
			}
		}
	}
	/*check OZ
		LOG(5,<<"sim matrix"<<endl);
		for (i=0;i<symmetric_q.size();++i) {
			for (j=0; j<symmetric_q.size(); ++j) {
				//LOG(5,<<symmetric_q[i][j]<<" ");
				LOG(5,<< setprecision(3) <<  setw(5) << symmetric_q[i][j]<<'\t');
				
			}
			LOG(5,<<endl);
			} */

}

void q2pt::calc_left_and_right_eig_of_pam(
		VVdouble &left_eig_of_pam,
		VVdouble &right_eig_of_pam,
		const VVdouble &v,
		const Vdouble& freq) {
//----------------------------------------------------------------------------------
//input:	left_eig_of_pam, right_eig_of_pam they will be the eigenvectors of pam1;	
//			freq is the vector of amino acid frequencies of the model.
//			v is the eigen vector matrix of the symetrical matrix
//output:	non
//doing:	now [SYM]  = [SqrtFreq] * [pam1] * inv([SqrtFreq])
//			so [pam1] = inv([SqrtFreq]) * [SYM] * [SqrtFreq]
//			SYM		  = [V] * [D] * transp([V])
//			hence [pam1] = {inv([SqrtFreq]) * [V]} * [D] * {transp([V]) * [SqrtFreq]}
//			{inv([SqrtFreq]) * [V]} is left_eig_of_pam, and the above one ^ is right.
//----------------------------------------------------------------------------------
	int i,j;
	for (i=0;i<v.size();++i) {
		for (j=0;j<v.size();++j)
		{
			if ((freq[i] != 0.0) &&(freq[j] != 0.0)) {
				left_eig_of_pam[i][j] =  (1/sqrt(freq[i]))* v[i][j];
				right_eig_of_pam[i][j]= sqrt(freq[j]) * v[j][i];
			}
		}
	}

//	LOG(5,<<"left_eig_of_pam"<<endl);
//	for (i=0;i<4;++i) {
//		for (j=0; j<4; ++j) {
//			LOG(5,<<left_eig_of_pam[i][j]<<" ");
//			LOG(5,<<pam1[i][i]<<" ");
//		}
//		LOG(5,<<endl);
//	}
//
//	LOG(5,<<"right eig_of_pam"<<endl);
//	for (i=0;i<4;++i) {
//		for (j=0; j<4; ++j) {
//			LOG(5,<<right_eig_of_pam[i][j]<<" ");
//			LOG(5,<<pam1[i][i]<<" ");
//		}
//		LOG(5,<<endl);
//	}
//
//	LOG(5,<<"press anykey"<<endl);
//	char lll;
//	cin>>lll;


}

VVdouble get1PamFromCountMatrix(const vector<MDOUBLE>& freq,
		   const VVdouble & sub_matrix){
//----------------------------------------------------------------------------------
//input:		pam1 : a pointer to the matrix where pam1 will be.
//				sub_matrix: the substitution matrix
//				freq vector: the amino acid's frequenceis.
//output:		non
//doing:		fill in 1 pam from sub matrix and freq vector
//calculation:  sub_matrix[a][b] is the substitution matrix, between a and b
//				(sub_matrix[a][b]=sub_matrix[b][a])
//				we use f[a][b] insted of sub_matrix[a][b] to be the same as the book
//(reference)	"introduction to computational molecular biology by setubal and meidanis pg 80;
//				let f[a] be sigma f[a][b] on all b (we made f[a][a] = 0;)
//				i.e. f[a] is the number of mutation from a observed
//				let f be sigma f[a] on all a; (=the total mutations*2)
//				now, the mutaibility of a is defined as 
//
//				(1)	m[a] = f[a] / (100*f*freq[a])
//
//				100*f is a scaling factor for 1 pam.
//				then pam1[a][b] will be pr(a->b/a changed) * pr(a changed)
//				
//				(2) pam1[a][b] = (f[a][b]/f[a])*m[a]
//
//				(3) f[a][a] = 1-m[a] (easy to show)
//
//				notice that sigma 1pam[a][b] over all b is 1 and that
//				sigma freq[a]*1pam[a][a] over all a is 0.99
//----------------------------------------------------------------------------------
	const int _alphabetSize=sub_matrix.size();
	VVdouble pam1;
	pam1.resize(_alphabetSize);
	for (int z=0; z < _alphabetSize; ++z) {
		pam1[z].resize(_alphabetSize,0);
	}

	int i,j;//indices
	MDOUBLE total=0;			// i.e.f in the above explanation
	for (i=0;i<_alphabetSize;++i) {
		for (j=0; j<_alphabetSize; ++j){
			total+=sub_matrix[i][j];
		}
	}
	
	MDOUBLE tmsum;
	for (i=0;i<_alphabetSize;++i) {
		tmsum = 0.0;
		for (j=i+1; j<_alphabetSize; ++j){
			if ((freq[i] == 0.0) || (freq[j] == 0.0)) {
				pam1[i][j] = 0.0;pam1[j][i] = 0.0;
			} else {
				pam1[i][j] = sub_matrix[i][j]/(100.0*total*freq[i]);
				pam1[j][i] = sub_matrix[i][j]/(100.0*total*freq[j]);
			}
		}
	}

	for (i=0;i<_alphabetSize;++i) {
		tmsum = 0.0;
		for (j=0;j<_alphabetSize;++j) {
			if (j!=i) tmsum += pam1[i][j];
		}

		if (freq[i] != 0.0)  {
			pam1[i][i]=1.0-tmsum;
		}
	}

#ifdef VERBOS
	LOG(5,<<" priting the 4*4 top-left corner of the 1pam matrix * 10^6 "<<endl);
	for (int a=0; a < 4; ++a) {
		for (int b=0; b < 4; ++b) {
			LOG(5,<<pam1[a][b]*1000000.0<<"   ");
		}
		LOG(5,<<endl);
	}
#endif
	return pam1;

}

