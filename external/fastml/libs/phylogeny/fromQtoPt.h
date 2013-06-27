// $Id: fromQtoPt.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___FROM_Q_TO_PT
#define ___FROM_Q_TO_PT

#include "replacementModel.h"
#include <cmath>
#include <iomanip>

int MyJacobi(VVdouble &Insym, VVdouble &RightEigenV, Vdouble &EigenValues);// num rec

VVdouble get1PamFromCountMatrix(const vector<MDOUBLE>& freq,
		const VVdouble & sub_matrix);

class q2pt : public replacementModel {
public:
	void fillFromRateMatrix(const vector<MDOUBLE>& freq,
					const VVdouble & qMatrix);
	void fillFrom1PAMMatrix(const vector<MDOUBLE>& freq,
					const VVdouble & onePam);

	
	explicit q2pt(): err_allow_for_pijt_function(1e-4){} 

	// @@@@ I'm not sure why I had to implement this operator=, but it doesn't work without it.
	q2pt& operator=(const q2pt &other) {
		_freq = other._freq;
		_leftEigen = other._leftEigen;
		_rightEigen = other._rightEigen;
		_eigenVector = other._eigenVector;
		return (*this);
	}

	virtual replacementModel* clone() const { return new q2pt(*this); }
//	virtual nucJC* clone() const { return new nucJC(*this); } // see note down:

	const int alphabetSize() const {return _freq.size();}


	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE freq(const int i) const {return _freq[i];};
	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE err_allow_for_pijt_function; //1e-4

	VVdouble getLeftEigen() const {return _leftEigen;} ;
	VVdouble getRightEigen() const {return _rightEigen;};
	Vdouble getEigenVec() const {return _eigenVector;};

private:
	Vdouble _freq;
	VVdouble _leftEigen;
	VVdouble _rightEigen;
	Vdouble _eigenVector;
	bool currectFloatingPointProblems(MDOUBLE& sum) const;

public: // to become private:
	void calc_symmetric_q(const VVdouble &q_matrix,VVdouble &symmetric_q,const Vdouble & freq);
	void calc_left_and_right_eig_of_pam(
		VVdouble &left_eig_of_pam,
		VVdouble &right_eig_of_pam,
		const VVdouble &v,
		const Vdouble& freq);
};

#endif

