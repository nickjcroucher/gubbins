#include "matrixUtils.h"
#include "errorMsg.h"
#include <cmath>
#include <string>
#include <ctype.h>
#include <cctype>
#include <cstdlib>

Vdouble getDiagonalFromMatrix(VVdouble &mat){
	Vdouble diagonal;
	for (int i=0; i<mat.size(); i++)
		diagonal.push_back(mat[i][i]);
	return diagonal;
}

Vdouble getSubDiagonalFromMatrix(VVdouble &mat){
	Vdouble diagonal;
	for (int i=0; i<mat.size()-1; i++)
		diagonal.push_back(mat[i+1][i]);
	return diagonal;
}






void readMatrixFromFile(VVdouble &mat,string fileName){
	ifstream in(fileName.c_str());
	if (!in){
		string err="in function readMatrixFromFile, empty file or non-existant:";
		err+=fileName;
		errorMsg::reportError(err); 
	}
	int i=0;
	mat.resize(1);
	while (!in.eof()) {
		string row;
		int k=0;
		getline(in,row,'\n');
		while (k<row.size()){
			string value;
			while (row[k]!=' ' && k<row.size()){
				value+=row[k];
				k++;
			}
		k++;
		mat[i].push_back(atof(value.c_str()));
		//j++;
		//mat.resize(j);
		}
	if (!in.eof())
		mat.resize(++i+1);
	} 
	in.close();
}


void printMatrix(const VVdouble &mat, ostream &out) {
	int num=mat.size();
	for (int row=0; row<num; row++) {
		for (int position=0; position<mat[row].size(); position++) {
			out << mat[row][position] << '\t';
		}
		out << endl ;
	}
	out << endl ;

}

void printMatrix(const VVint &mat, ostream &out) {
	int num=mat.size();
	for (int row=0; row<num; row++) {
		for (int position=0; position<mat[row].size(); position++) {
			out << mat[row][position] << '\t';
		}
		out << endl ;
	}
	out << endl ;
	out<<"---------------------------------------------"<<endl;

}



VVdouble transpose(const VVdouble &mat){
	VVdouble matT;
	int n=mat.size();
	resizeMatrix(matT,n,n);
	for (int i=0; i<n;i++){
		for (int j=0; j<n;j++) {
			matT[i][j]=mat[j][i];
		}
	}
	return matT;
}




VVdouble subtract(const VVdouble &mat1,const VVdouble &mat2){
	VVdouble newMat=add(mat1,reverseSign(mat2));
	return newMat;
}

VVdouble reverseSign(const VVdouble &mat1){
	VVdouble newMat(mat1.size());
	for (int i=0;i<mat1.size();i++){
		newMat[i].resize(mat1[i].size());
		for (int j=0;j<mat1.size();j++){
			newMat[i][j]=-mat1[i][j];
		}
	}
	return newMat;

}


void findMaxInVector(const Vdouble &vec, MDOUBLE &maxValue, int &argmax){
	MDOUBLE tempMax=VERYSMALL;
	int tempArgMax=0;
	for (int i=0; i<vec.size(); i++){
		if (vec[i]>tempMax){
			tempMax=vec[i];
			tempArgMax=i;
		}
	}
	maxValue=tempMax;
	argmax=tempArgMax;
}

void findMinInVector(const Vdouble &vec, MDOUBLE &minValue, int &argmin) {
	Vdouble minusCopy(vec.size());
	for (int i=0; i<vec.size(); i++){
		minusCopy[i] = -vec[i];
	}
	findMaxInVector(minusCopy, minValue, argmin);
	minValue = -minValue;
}

MDOUBLE averageElementInVector(const Vdouble &vec) {
	MDOUBLE sum=0.0;
	for (int i=0; i<vec.size(); i++){
		sum+=vec[i];
	}
	return sum/vec.size();
}

void appendBinaryVectors(Vint &vec1, const Vint &vec2){
    for (int i=0; i < vec2.size(); i++) 
        if (vec2[i]==1)
			vec1[i]=1;
}

void appendVectors(Vint &vec1, const Vint &vec2) {
	for (int i=0; i<vec2.size();i++)
		vec1.push_back(vec2[i]);
}

Vint complementBinaryVec(Vint&bufferVec) {
	for (int i=0; i<bufferVec.size(); i++) 
		bufferVec[i]=abs(bufferVec[i]-1);
	return bufferVec;
}


//reads a vertical vector of float numbers(separated by \n)
void readDoubleVecFromFile(Vdouble &vec,string fileName){
	ifstream in(fileName.c_str());
	if (!in){
		string err="in function readDoubleVecFromFile, empty file or non-existant:";
		err+=fileName;
		errorMsg::reportError(err); 
	}
	string row;
	while (!in.eof()){
		getline(in,row,'\n');
		//if (isalnum(*(row.c_str())) || (row[0]=="."))
		if (isspace(*(row.c_str())) || row=="") continue;
			vec.push_back(atof(row.c_str()));
	}

	in.close();
}

void normalize(Vdouble &vec){
	MDOUBLE sum=0.0;
	MDOUBLE squareSum=0.0;
	int N=vec.size();
	int i=0;
	for (i=0;i<N;i++) sum+=vec[i];
	for (i=0;i<N;i++) squareSum+=(vec[i]*vec[i]);
	MDOUBLE avg=sum/N;
	MDOUBLE sqrAvg=squareSum/N;
	MDOUBLE stdDev=sqrt(sqrAvg-avg*avg);
	for (i=0;i<N;i++) vec[i]=(vec[i]-avg)/stdDev;

}

void scaleByAverage(Vdouble &vec){
	MDOUBLE sum=0.0;
	MDOUBLE squareSum=0.0;
	int N=vec.size();
	int i=0;
	for (i=0;i<N;i++) sum+=vec[i];
	for (i=0;i<N;i++) squareSum+=(vec[i]*vec[i]);
	MDOUBLE avg=sum/N;
	for (i=0;i<N;i++) vec[i]=(vec[i])/avg;
}

Vdouble solveLinearEquations(VVdouble A,Vdouble b){
//	VVdouble Acopy=A; //creating a copy, since ludcmp&lubksb destroy the input
//	Vdouble bcopy=b;
	MDOUBLE d; //required for ludcmp; irrelevant for us.
	Vdouble indx; //required for ludcmp; irrelevant for us.
	ludcmp(A,indx,d); //decomposes A into product of diagonal matrices
	lubksb(A,indx,b); //solves
	return b;
}


void ludcmp(VVdouble &a, Vdouble &indx, MDOUBLE &d)
{
	const MDOUBLE TINY=1.0e-20;
	int i,imax=0,j,k;
	MDOUBLE big,dum,sum,temp;

	int n=a.size();
	Vdouble vv(n);
	indx.resize(n);//my addition
	d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) errorMsg::reportError("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ((dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			d = -d;
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
}



void lubksb(VVdouble &a, Vdouble &indx, Vdouble &b)
{
	int i,ii=0,ip,j;
	MDOUBLE sum;

	int n=a.size();
	for (i=0;i<n;i++) {
		ip=(int)(indx[i]);
		sum=b[ip];
		b[ip]=b[i];
		if (ii != 0)
			for (j=ii-1;j<i;j++) sum -= a[i][j]*b[j];
		else if (sum != 0.0)
			ii=i+1;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

//get the first norm sum{abs(Mij)}
MDOUBLE getMatrixNorm(const VVdouble &mat) {
	MDOUBLE res(0.0);
	for (int i=0; i<mat.size(); i++){
		for (int j=0; j<mat[i].size();j++){
			res += fabs(mat[i][j]);
		}
	}
	return res;
}

/********************************************************************************************
*********************************************************************************************/
void resize_VVVV(int dim1, int dim2, int dim3, int dim4,  VVVVdouble& vetor){
	
	vetor.resize(dim1);
	for (int posNum=0;posNum<vetor.size();++posNum){
		vetor[posNum].resize(dim2);
		for (int n=0;n<vetor[posNum].size();++n){
			resizeMatrix(vetor[posNum][n],dim3,dim4);
		}
	}
}
/********************************************************************************************
*********************************************************************************************/
void resize_VVV(int dim1, int dim2, int dim3, VVVdouble& vetor){	
	vetor.resize(dim1);
	for (int n=0;n<vetor.size();++n){
		resizeMatrix(vetor[n],dim2,dim3);
	}
}



