#ifndef ___MATRIX_UTIL_H
#define ___MATRIX_UTIL_H

#include "definitions.h"
#include "logFile.h"
#include "errorMsg.h"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

class sequenceContainer;
using namespace std;



void printMatrix(const VVdouble &mat, ostream &out);
void printMatrix(const VVint &mat, ostream &out) ;
 
void readMatrixFromFile(VVdouble &mat,string fileName);

Vdouble getDiagonalFromMatrix(VVdouble &mat);
Vdouble getSubDiagonalFromMatrix(VVdouble &mat);

//get the first norm sum{abs(Mij)}
MDOUBLE getMatrixNorm(const VVdouble &mat);


template<typename _T>
void unitMatrix(vector<vector< _T> > &m, int n){
	resizeMatrix(m,n,n);
	for (int i=0; i<n; i++){
		for (int j=0; j<n;j++){
			if (i==j) m[i][j]=1;
			else m[i][j]=0;
		}
	}
}

template<typename _T>
void zeroMatrix(vector<vector< _T> > &m){
	for (int i=0; i < m.size(); i++)
		for (int j=0; j<m[i].size();j++)
			m[i][j]=0;
}

template<typename _T>
void oneMatrix(vector<vector< _T> > &m){
	for (int i=0; i < m.size(); i++)
		for (int j=0; j<m[i].size();j++)
			m[i][j]=1;
}


//assumes that #columns in mat1=#rows in mat2
template<typename _T>
vector<vector< _T> > multiplyMatrixes(vector<vector< _T> > &mat1, vector<vector< _T> > &mat2){
	vector<vector< _T> >  mat;
	if ((mat1.size()==0) || (mat2.size() ==0))
		errorMsg::reportError("Error in multiplyMatrixes, one of the matrices inputted is of size 0");;
	int numColumns=mat1[0].size();
	int numRows = mat2.size();
	resizeMatrix(mat,numColumns,numRows);
	for (int i=0; i<numColumns; i++){
		for (int j=0; j<numRows;j++){
			for (int k=0;k<numColumns;k++){
				mat[i][j]+=mat1[i][k]*mat2[k][j];
			}
		}
	}
	return mat;
}

template<typename _T>
vector<vector< _T> > multiplyMatrixByScalar(const vector<vector< _T> > &mat, MDOUBLE scalar) {
	vector<vector< _T> > mat_copy = mat;
	for (int i=0; i<mat.size(); i++){
		for (int j=0; j<mat[i].size();j++){
			mat_copy[i][j]*=scalar;
		}
	}
	return mat_copy;
}

template<typename _T>
void resizeMatrix(vector<vector< _T> > &mat, int rows, int columns){
	mat.resize(rows);
	for (int i=0; i<rows;i++){
		mat[i].resize(columns);
		for (int j=0;j<columns;j++){ // initializing all values as zero
			mat[i][j] = 0.0;
		}
	}
}

template<typename _T>
vector<vector< _T> > add(const vector<vector< _T> >  &mat1,const vector<vector< _T> >  &mat2){
	if (mat1.size()!=mat2.size()) errorMsg::reportError("different sized matrices in matrixUtils::add");
	vector<vector< _T> > newMat(mat1.size());
	for (int i=0;i<mat1.size();i++){
		if (mat1[i].size()!=mat2[i].size()) errorMsg::reportError("different sized matrices in matrixUtils::add");
		newMat[i].resize(mat1[i].size());
		for (int j=0;j<mat1.size();j++){
			newMat[i][j]=mat1[i][j]+mat2[i][j];
		}
	}
	return newMat;
}

template<typename _T> 
void printVec(vector< _T> &vec,ostream &out=cout,bool printVertical=true) {
	for (int i=0; i<vec.size();i++){
		out<< vec[i];
		out<<(printVertical?"\n":" ");
	}
	out<<endl;
}



VVdouble transpose(const VVdouble &mat);
VVdouble subtract(const VVdouble &mat1,const VVdouble &mat2);
VVdouble reverseSign(const VVdouble &mat1);

void findMaxInVector(const Vdouble &vec, MDOUBLE &maxValue, int &argmax)  ;
void findMinInVector(const Vdouble &vec, MDOUBLE &minValue, int &argmin)  ;
MDOUBLE averageElementInVector(const Vdouble &vec)  ;
void appendBinaryVectors(vector <int> &vec1, const vector <int> &vec2);
void appendVectors(Vint &vec1, const Vint &vec2);
Vint complementBinaryVec(vector <int>&bufferVec) ; // returns complementary binary vector
void readDoubleVecFromFile(Vdouble &vec,string fileName); //reads a vertical vector (separated by \n)

void normalize(Vdouble &vec);
void scaleByAverage(Vdouble &vec);


//solve nxn linear equations of the form Ax=b; return x;
Vdouble solveLinearEquations(VVdouble A,Vdouble b);
// functions from numerical recipes that solve nxn linear equations
void lubksb(VVdouble &a, Vdouble &indx, Vdouble &b);
void ludcmp(VVdouble &a, Vdouble &indx, MDOUBLE &d);

void resize_VVVV(int dim1, int dim2, int dim3, int dim4,  VVVVdouble& vetor);
void resize_VVV(int dim1, int dim2, int dim3, VVVdouble& vetor);




#endif
