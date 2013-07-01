#include "wYangModel.h"
#include "codon.h"
#include "readDatMatrix.h" // for the normalizeQ function.

wYangModel::wYangModel(const MDOUBLE inW, const MDOUBLE inK,bool globalW,  codon * coAlph):
	_w(inW),_k(inK),_globalW(globalW),_coAlpha(NULL){
	_coAlpha = (codon*)(coAlph->clone());
	codonUtility::initSubMatrices(*_coAlpha);
	homogenousFreq();
	_Q.resize(alphabetSize());
	for (int z=0; z < _Q.size();++z) _Q[z].resize(alphabetSize(),0.0);
	updateQ();
}

wYangModel::wYangModel(const MDOUBLE inW, const MDOUBLE inK, const Vdouble& freq,bool globalW, codon * coAlph):
	_w(inW),_k(inK),_globalW(globalW),_freq(freq),_coAlpha(NULL){
	_coAlpha = (codon*)(coAlph->clone());
	_Q.resize(alphabetSize());
	codonUtility::initSubMatrices(*_coAlpha);
	for (int z=0; z < _Q.size();++z) _Q[z].resize(alphabetSize(),0.0);
	updateQ();
}


wYangModel& wYangModel::operator=(const wYangModel &other) {
	_w = other._w;
	_k = other._k;
	_q2pt = other._q2pt;
	_Q = other._Q;
	_globalW = other._globalW;
	_freq = other._freq;
	if (_coAlpha) delete _coAlpha;
	if (other._coAlpha) 
		_coAlpha = (codon*)(other._coAlpha->clone());
	else
		_coAlpha = NULL;
	return *this;

}



void wYangModel::updateQ() {
	int i,j;
	MDOUBLE sum=0.0;
	for (i=0; i < _Q.size();++i) {
		for (j=i+1; j < _Q.size();++j) {
			MDOUBLE val;
			if (codonUtility::codonReplacement(i,j) == codonUtility::non_synonymous) {
				if (codonUtility::codonDiff(i,j) == codonUtility::tr) val = _k*_w;	
				else if (codonUtility::codonDiff(i,j) == codonUtility::tv) val = _w;
				else val = 0;//more than one substitution.
			}
			else {//synonymous
				if (codonUtility::codonDiff(i,j) == codonUtility::tr) val = _k;	
				else if (codonUtility::codonDiff(i,j) == codonUtility::tv) val = 1;
				else val = 0;//more than one substitution.
			}
			_Q[i][j] = val * _freq[j];
			_Q[j][i] = val * _freq[i];
		}
		_Q[i][i] = 0.0;  //temporary value
	}
	// filling the diagonal
	for (i=0; i < _Q.size(); ++i){
		sum = 0.0;
		for (j=0; j < _Q.size(); ++j) {
			sum += _Q[i][j];
		}
		_Q[i][i] = -sum;
	}
	if (_globalW == true) // w is not distributed, only one Q matrix
		normalizeQ(_Q,_freq);

	_q2pt.fillFromRateMatrix(_freq,_Q); 
}


void wYangModel::norm(MDOUBLE scale){
	for (int i=0; i < _Q.size(); ++i) {
		for (int j=0; j < _Q.size(); ++j) {
			_Q[i][j] *=scale; 
			
		}
	}
	_q2pt.fillFromRateMatrix(_freq,_Q);
}


MDOUBLE wYangModel::sumPijQij(){
	MDOUBLE sum=0.0;
	for (int i=0; i < _Q.size(); ++i) {
		sum -= (_Q[i][i])*_freq[i];
	}
	return sum;
}
