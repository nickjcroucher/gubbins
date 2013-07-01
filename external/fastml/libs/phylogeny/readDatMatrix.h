// $Id: readDatMatrix.h 5805 2009-01-20 09:19:26Z adido $

#ifndef ___READ_DAT_MATRIX
#define ___READ_DAT_MATRIX

#include "definitions.h"
#include <string>
#include <iostream>
#include <fstream>
#include "datMatrixHolder.h"

using namespace std;

void normalizeQ(VVdouble& q, const Vdouble& freq);

void readDatMatrixFromFile(const string & matrixFileName,
						   VVdouble & subMatrix,
						   Vdouble & freq);
void readDatMatrixFromString(const string & matrixFileString,
			     VVdouble & subMatrix,
			     Vdouble & freq, int alphaSize = 20);

VVdouble fromWagSandFreqToQ(const VVdouble & s,const Vdouble& freq);

#include "replacementModel.h"
#include "definitions.h"
#include "errorMsg.h"

class pupAll : public replacementModel {
public:
	// get matrix from file:
	explicit pupAll(const string& matrixFileString) : err_allow_for_pijt_function(1e-4) {fillMatricesFromFile(matrixFileString);}
	explicit pupAll(const string& matrixFileString, const vector<MDOUBLE>& freq) : err_allow_for_pijt_function(1e-4) {fillMatricesFromFile(matrixFileString,freq);}

	// get matrix from within the .exe
	explicit pupAll(const datMatrixString& matrixFileString,int alphaSize = 20) : err_allow_for_pijt_function(1e-4) {fillMatrices(matrixFileString.Val,alphaSize); }
	explicit pupAll(const datMatrixString& matrixFileString, const vector<MDOUBLE>& freq) : err_allow_for_pijt_function(1e-4) {fillMatrices(matrixFileString.Val,freq);}


	const int alphabetSize() const {return _freq.size();}//20 or 61
	const MDOUBLE err_allow_for_pijt_function; //1e-4
	virtual replacementModel* clone() const { return new pupAll(*this); }

	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE t) const;
	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE t) const;
	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE t) const;
	const MDOUBLE freq(const int i) const {return _freq[i];}

	const MDOUBLE Pij_tAlpha    (const int i,const int j, const MDOUBLE t, const MDOUBLE alpha) const;
	const MDOUBLE Pij_tAlpha_dt (const int i,const int j, const MDOUBLE t, const MDOUBLE alpha) const;
	const MDOUBLE Pij_tAlpha_dt2(const int i,const int j, const MDOUBLE t, const MDOUBLE alpha) const;

private:
	void fillMatrices(const string & matrixName,const vector<MDOUBLE>& freq);
	void fillMatrices(const string & matrixName,int alphaSize);
	void fillMatricesFromFile(const string & dataFileString,const vector<MDOUBLE>& freq);
	void fillMatricesFromFile(const string & dataFileString);


	bool currectFloatingPointProblems(MDOUBLE& sum) const;

	VVdouble _leftEigen;
	VVdouble _rightEigen;
	Vdouble _eigenVector;
	Vdouble _freq;
};

#endif
