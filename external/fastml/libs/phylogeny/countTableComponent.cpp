// $Id: countTableComponent.cpp 962 2006-11-07 15:13:34Z privmane $

// version 1.00
// last modified 3 Nov 2002

#include "countTableComponent.h"
#include "logFile.h"

void countTableComponentHom::zero() {
	for (int i=0; i < _countValues.size() ;++i) {
		for (int j=0; j < _countValues[0].size() ;++j) {
				_countValues[i][j] = 0;
		}
	}
}

void countTableComponentHom::countTableComponentAllocatePlace(
		const int alphabetSize) {
	int i;
	_countValues.resize(alphabetSize);
	for (i=0; i < alphabetSize;++i) _countValues[i].resize(alphabetSize);
}

void countTableComponentHom::printTable(ostream& out) const {
	MDOUBLE sumCheck = 0.0;
	for (int i=0; i < _countValues.size();++i) {
		for (int k=0; k <  _countValues.size();++k) {
			out<<"counts["<<i<<"]["<<k<<"]"<<_countValues[i][k];
			sumCheck += _countValues[i][k];
			out<<endl;
		}
	}
	out<<"sum is: "<<sumCheck<<endl;
}

