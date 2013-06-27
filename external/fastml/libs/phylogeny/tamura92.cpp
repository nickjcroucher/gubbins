// $Id: tamura92.cpp 962 2006-11-07 15:13:34Z privmane $

#include "tamura92.h"
#include "errorMsg.h"

// This implementation was copied from the Bio++ Phyl library (by Julien Dutheil) - file T92.cpp

tamura92::tamura92(const MDOUBLE theta,
				   const MDOUBLE TrTv) 
	: _theta(theta), _TrTv(TrTv) {

	_freq.resize(4);
	changeTheta(theta);
}

void tamura92::changeTheta(const MDOUBLE theta) {
	_theta = theta;
	_freq[0] = _freq[3] = (1.0 - theta) / 2.0;
	_freq[1] = _freq[2] = theta / 2.0;
}

const MDOUBLE tamura92::Pij_t(const int i, const int j, const MDOUBLE t) const {
	double k = (_TrTv + 1.0) / 2.0;
	double r  = 2.0 / (1.0 + 2.0 * _theta * _TrTv - 2.0 * _theta * _theta * _TrTv);
	double l = r * t;
	double exp1 = exp(-l);
	double exp2 = exp(-k * l);
	
	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return _freq[0] * (1.0 + exp1) + _theta * exp2; //A
				case 1 : return _freq[1] * (1.0 - exp1);                 //C
				case 2 : return _freq[2] * (1.0 + exp1) - _theta * exp2; //G
				case 3 : return _freq[3] * (1.0 - exp1);                 //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return _freq[0] * (1.0 - exp1);                        //A
				case 1 : return _freq[1] * (1.0 + exp1) + (1. - _theta) * exp2; //C
				case 2 : return _freq[2] * (1.0 - exp1);                        //G
				case 3 : return _freq[3] * (1.0 + exp1) - (1. - _theta) * exp2; //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return _freq[0] * (1.0 + exp1) - (1. - _theta) * exp2; //A
				case 1 : return _freq[1] * (1.0 - exp1);                        //C
				case 2 : return _freq[2] * (1.0 + exp1) + (1. - _theta) * exp2; //G
				case 3 : return _freq[3] * (1.0 - exp1);                        //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return _freq[0] * (1.0 - exp1);                 //A
				case 1 : return _freq[1] * (1.0 + exp1) - _theta * exp2; //C
				case 2 : return _freq[2] * (1.0 - exp1);                 //G
				case 3 : return _freq[3] * (1.0 + exp1) + _theta * exp2; //T, U
			}
		}
	}
	return -1;
}

const MDOUBLE tamura92::dPij_dt(const int i,const int j, const MDOUBLE t) const {
	double k = (_TrTv + 1.0) / 2.0;
	double r  = 2.0 / (1.0 + 2.0 * _theta * _TrTv - 2.0 * _theta * _theta * _TrTv);
	double l = r * t;
	double exp1 = exp(-l);
	double exp2 = exp(-k * l);

	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return r * (_freq[0] * - exp1 + _theta * -k * exp2); //A
				case 1 : return r * (_freq[1] *   exp1);                      //C
				case 2 : return r * (_freq[2] * - exp1 - _theta * -k * exp2); //G
				case 3 : return r * (_freq[3] *   exp1);                      //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return r * (_freq[0] *   exp1);                               //A
				case 1 : return r * (_freq[1] * - exp1 + (1.0 - _theta) * -k * exp2); //C
				case 2 : return r * (_freq[2] *   exp1);                               //G
				case 3 : return r * (_freq[3] * - exp1 - (1.0 - _theta) * -k * exp2); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return r * (_freq[0] * - exp1 - (1.0 - _theta) * -k * exp2); //A
				case 1 : return r * (_freq[1] *   exp1);                               //C
				case 2 : return r * (_freq[2] * - exp1 + (1.0 - _theta) * -k * exp2); //G
				case 3 : return r * (_freq[3] *   exp1);                               //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return r * (_freq[0] *   exp1);                      //A
				case 1 : return r * (_freq[1] * - exp1 - _theta * -k * exp2); //C
				case 2 : return r * (_freq[2] *   exp1);                      //G
				case 3 : return r * (_freq[3] * - exp1 + _theta * -k * exp2); //T, U
			}
		}
	}
	return -1;
}

const MDOUBLE tamura92::d2Pij_dt2(const int i,const int j, const MDOUBLE t) const {
	double k = (_TrTv + 1.0) / 2.;
	double k2 = k * k;
	double r  = 2.0 / (1.0 + 2.0 * _theta * _TrTv - 2.0 * _theta * _theta * _TrTv);
	double l = r * t;
	double r2 = r * r;
	double exp1 = exp(-l);
	double exp2 = exp(-k * l);

	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return r2 * (_freq[0] *   exp1 + _theta * k2 * exp2); //A
				case 1 : return r2 * (_freq[1] * - exp1);                      //C
				case 2 : return r2 * (_freq[2] *   exp1 - _theta * k2 * exp2); //G
				case 3 : return r2 * (_freq[3] * - exp1);                      //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return r2 * (_freq[0] * - exp1);                               //A
				case 1 : return r2 * (_freq[1] *   exp1 + (1.0 - _theta) * k2 * exp2); //C
				case 2 : return r2 * (_freq[2] * - exp1);                               //G
				case 3 : return r2 * (_freq[3] *   exp1 - (1.0 - _theta) * k2 * exp2); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return r2 * (_freq[0] *   exp1 - (1.0 - _theta) * k2 * exp2); //A
				case 1 : return r2 * (_freq[1] * - exp1);                               //C
				case 2 : return r2 * (_freq[2] *   exp1 + (1.0 - _theta) * k2 * exp2); //G
				case 3 : return r2 * (_freq[3] * - exp1);                               //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return r2 * (_freq[0] * - exp1);                      //A
				case 1 : return r2 * (_freq[1] *   exp1 - _theta * k2 * exp2); //C
				case 2 : return r2 * (_freq[2] * - exp1);                      //G
				case 3 : return r2 * (_freq[3] *   exp1 + _theta * k2 * exp2); //T, U
			}
		}
	}
	return -1;
}

