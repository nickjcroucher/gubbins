#ifdef DOUBLEREP
#include "doubleRep.h"
#include <cmath>



doubleRepMantisa::doubleRepMantisa(MDOUBLE mantissa, int expon){
	_mantissa=mantissa;
	_expon=expon;
	fixParams();
}


doubleRepMantisa::doubleRepMantisa(MDOUBLE a){
	int answerExp;
	MDOUBLE answerMantissa=frexp(a,&answerExp);
	_mantissa=answerMantissa;
	_expon=answerExp;
}



doubleRepMantisa::doubleRepMantisa(const doubleRepMantisa& other): _mantissa(other._mantissa), _expon(other._expon) {
}


//make sure 0.5<=mantissa<1, as a matter of convention
void doubleRepMantisa::fixParams(){
	while (_mantissa>=1){
		_expon++;
		_mantissa/=2.0;
	}
	while ((_mantissa<0.5) && (_mantissa>0)){
		_expon--;
		_mantissa*=2.0;
	}
	while (_mantissa<=-1){
		_expon++;
		_mantissa/=2.0;
	}
	while ((_mantissa>-0.5) && (_mantissa<0)){
		_expon--;
		_mantissa*=2.0;
	}
}

MDOUBLE convert(const doubleRepMantisa& a){
	MDOUBLE aFullRep= ldexp(a._mantissa,a._expon);
	return aFullRep;
}

//switches from base 2 to base e
const MDOUBLE doubleRepMantisa::d_log() const{
  static const MDOUBLE log2(log(2.0));
	return log(_mantissa)+log2*_expon;
}


ostream& operator<<(ostream &out, const doubleRepMantisa& a){
	a.output(out);
	//  a.output0x(out);
//	out<<a._mantissa<<string(" * 2^")<<a._expon;
//	out<<a._mantissa<<" * 2^"<<a._expon;
    return out;
}

istream& operator>>(istream &in, doubleRepMantisa& a) {
  MDOUBLE num;
  in >> num;
  a = num;
  return in;
}
#endif
