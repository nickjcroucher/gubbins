// $Id: goldmanYangModel.cpp 962 2006-11-07 15:13:34Z privmane $

#include "goldmanYangModel.h"
#include "codon.h"
#include "readDatMatrix.h" // for the normalizeQ function.


goldmanYangModel::goldmanYangModel(const MDOUBLE inV, const MDOUBLE inK,codon & inCodonAlph, const bool globalV):
		_v(inV),_k(inK),_globalV(_globalV),_codonAlph(inCodonAlph){
		homogenousFreq();
		_Q.resize(_codonAlph.size());
		for (int z=0; z < _Q.size();++z) _Q[z].resize(_codonAlph.size(),0);
		updateQ();
	
}


goldmanYangModel::goldmanYangModel(const MDOUBLE inV, const MDOUBLE inK, codon & inCodonAlph,const Vdouble& freq,const bool globalV):
		_freq(freq),_v(inV),_k(inK),_globalV(_globalV),_codonAlph(inCodonAlph){
		_Q.resize(_codonAlph.size());
		for (int z=0; z < _Q.size();++z) _Q[z].resize(_codonAlph.size(),0);
		updateQ();
}


void goldmanYangModel::updateQ() {
	
	// building q.
	int i,j;
	MDOUBLE sum=0.0;
	MDOUBLE epsilon=0.00000001;//0.00000000001;
	MDOUBLE factor = 1000.0; 
	for (i=0; i < _Q.size();++i) {
		sum=0;
		for (j=0; j < _Q.size();++j) {
			if (j==i) continue; //same codon
			if (codonUtility::codonDiff(i,j,_codonAlph) == codonUtility::tr) {
				_Q[i][j] = _k*exp(-(1/factor)*_gcd.getGranthamDistance(codonUtility::aaOf(i,_codonAlph),codonUtility::aaOf(j,_codonAlph))*_v);	
				if (_Q[i][j]<epsilon) _Q[i][j] = epsilon;
			}else if (codonUtility::codonDiff(i,j,_codonAlph) == codonUtility::tv) {
				_Q[i][j] = exp(-(1/factor)*_gcd.getGranthamDistance(codonUtility::aaOf(i,_codonAlph),codonUtility::aaOf(j,_codonAlph))*_v);
				if (_Q[i][j]<epsilon) _Q[i][j] = epsilon;
			}
			else _Q[i][j] = 0;//more than one substitution.
			
			_Q[i][j]*=_freq[j];
			sum += _Q[i][j];

		}
		_Q[i][i]=-sum;
	}
	

	// check:
/*	LOG(5,<<"\n\n\n ===================================== \n");
	int a1,a2;
	for (a1=0;a1<4;++a1){
		for (a2=0;a2<4;++a2){
			LOG(5,<<qMatrix[a1][a2]<<"\t");
		}
		LOG(5,<<endl);
	}
*/


	if (_globalV == true)
		normalizeQ(_Q,_freq);
	
	// check:
/*	LOG(5,<<"\n\n\n ===================================== \n");
	for (a1=0;a1<4;++a1){
		for ( a2=0;a2<4;++a2){
			LOG(5,<<qMatrix[a1][a2]<<"\t");
		}
		LOG(5,<<endl);
	}
*/

	
	// updating _q2Pt;
//	_Q = qMatrix;
	_q2pt.fillFromRateMatrix(_freq,_Q); 

	
	
}



// original with V and not 1/V 
/*
void goldmanYangModel::updateQ() {
	// building q.
	VVdouble qMatrix(_codonAlph.size());
	int i,j,z;
	MDOUBLE sum=0.0;
	for (z=0; z < qMatrix.size();++z) qMatrix[z].resize(_codonAlph.size(),0);
	for (i=0; i < qMatrix.size();++i) {
		sum=0;
		for (j=0; j < qMatrix.size();++j) {
			if (j==i) continue;
			if (codonUtility::codonDiff(i,j) == codonUtility::different) {
				qMatrix[i][j] =0;
			} else if (codonUtility::codonDiff(i,j) == codonUtility::transition) {
				qMatrix[i][j] =_k*exp(-_gcd.getGranthamDistance(codonUtility::aaOf(i),codonUtility::aaOf(j))/_v);
			} else if (codonUtility::codonDiff(i,j) == codonUtility::transversion) {
				qMatrix[i][j] = exp(-_gcd.getGranthamDistance(codonUtility::aaOf(i),codonUtility::aaOf(j))/_v);
			}
			qMatrix[i][j]*=_freq[j];
			sum += qMatrix[i][j];
		}
		qMatrix[i][i]=-sum;
	}
	// check:
	//LOG(5,<<"\n\n\n ===================================== \n");
	//int a1,a2;
	//for (a1=0;a1<4;++a1){
	//	for (a2=0;a2<4;++a2){
	//		LOG(5,<<qMatrix[a1][a2]<<"\t");
	//	}
	//	LOG(5,<<endl);
	//}

	if (_globalV == true)
		normalizeQ(qMatrix,_freq);

	//LOG(5,<<"\n\n\n ===================================== \n");
	//LOG(5,<<endl<<endl);
	//for (a1=0;a1<4;++a1){
	//	for (a2=0;a2<4;++a2){
	//		LOG(5,<<qMatrix[a1][a2]<<"\t");
	//	}
	//	LOG(5,<<endl);
	//}
	
	// updating _q2Pt;
	_Q = qMatrix;
	_q2pt.fillFromRateMatrix(_freq,qMatrix);
}


*/


