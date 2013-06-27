#ifndef __DOUBLE_REP_H
#define __DOUBLE_REP_H

#ifdef DOUBLEREP
#include "definitions.h"

#include <iostream>
#include <cmath>
using namespace std;

/* doubleRepMantisa: enables working with much larger or smaller numbers than normally possible
by the regular double representation
 * Representation of a double x as x=_mantissa*2^_expon
 Note: Base is 2!!
 */ 

class doubleRepMantisa{
public:
	
	doubleRepMantisa(){};
	explicit doubleRepMantisa(MDOUBLE mantissa, int expon);
	doubleRepMantisa(MDOUBLE a);
	doubleRepMantisa(const doubleRepMantisa& other);
	doubleRepMantisa* clone() {return new doubleRepMantisa(*this);}

	void output(ostream &out) const{ out<<_mantissa<<string(" * 2^")<<_expon;}
 //   void output0x(ostream &out) const{ double e0x=_expon*0.3010299956639; // log_10(2)
	//  int e=(int)(trunc(e0x))-1; 
	//  double m=_mantissa*pow(10,e0x-e); 
	//  out<<m;
	//  if (e<0)
	//	out<<"e"<<e;
	//  else
	//	out<<"e+"<<e;
	//}
	void outputn(ostream &out) { out<<_mantissa<<string(" * 2^")<<_expon<<endl;}
	
	friend MDOUBLE convert(const doubleRepMantisa& a);
	inline doubleRepMantisa& operator=(const doubleRepMantisa& a);
	inline doubleRepMantisa& operator+=(doubleRepMantisa a);
    inline doubleRepMantisa& operator++(); 
    inline doubleRepMantisa operator++(int);
    inline doubleRepMantisa& operator--(); 
    inline doubleRepMantisa operator--(int);
	friend inline doubleRepMantisa operator+(const doubleRepMantisa& a, const doubleRepMantisa& b);
	inline doubleRepMantisa& operator-=(const doubleRepMantisa& a);
	friend inline doubleRepMantisa operator-(const doubleRepMantisa& a, const doubleRepMantisa& b);
	inline doubleRepMantisa& operator*=(const doubleRepMantisa& a);
	friend inline doubleRepMantisa operator*(const doubleRepMantisa& a, const doubleRepMantisa& b);
	inline doubleRepMantisa& operator/=(const doubleRepMantisa& a);
	friend inline doubleRepMantisa operator/(const doubleRepMantisa& a, const doubleRepMantisa& b);

	friend inline bool operator==(const doubleRepMantisa& a, const doubleRepMantisa& b);
	friend inline bool operator!=(const doubleRepMantisa& a, const doubleRepMantisa& b);
	friend inline bool operator<(const doubleRepMantisa& a, const doubleRepMantisa& b);
	friend inline bool operator<=(const doubleRepMantisa& a, const doubleRepMantisa& b);
	friend inline bool operator>(const doubleRepMantisa& a, const doubleRepMantisa& b);
	friend inline bool operator>=(const doubleRepMantisa& a, const doubleRepMantisa& b);
	friend inline doubleRepMantisa abs(const doubleRepMantisa& d);

	
    const MDOUBLE d_log() const;
//	friend ostream& operator<<(ostream &out, const doubleRepMantisa& a);

	const MDOUBLE mantissa() const {return _mantissa;}
	const int expon() const {return _expon;}

private:
	void fixParams(); 


private:
	MDOUBLE _mantissa;
	int _expon;
};
	
inline doubleRepMantisa& doubleRepMantisa::operator=(const doubleRepMantisa& a){
	_mantissa=a.mantissa();
	_expon=a.expon();
	return *this;
}


inline doubleRepMantisa& doubleRepMantisa::operator++() {
  return (*this)+=1;
}

// matan:
inline doubleRepMantisa doubleRepMantisa::operator++(int) {
  doubleRepMantisa ans = *this;
  ++(*this);
  return ans;
}

// matan:
inline doubleRepMantisa& doubleRepMantisa::operator--() {
  return (*this)-=1;
}

// matan:
inline doubleRepMantisa doubleRepMantisa::operator--(int) {
  doubleRepMantisa ans = *this;
  --(*this);
  return ans;
}


// Original version by Adi Stern
inline doubleRepMantisa& doubleRepMantisa::operator+=(doubleRepMantisa a){
	//ensuring that (*this) is bigger than 'a' for sake of convenience
	if (a.expon()>_expon || ((a.expon()==_expon) && (a.mantissa()>_mantissa))){
		MDOUBLE tmpMant=0.0; int tmpExp=0;
		tmpMant=_mantissa;
		tmpExp=_expon;
		_mantissa=a.mantissa();
		a._mantissa=tmpMant;
		tmpExp=_expon;
		_expon=a.expon();
		a._expon=tmpExp;
	}
	if (a.mantissa()==0)
		return *this;
	if (_mantissa==0){
		_mantissa=a.mantissa();
		_expon=a.expon();
		return *this;
	}
	if (abs(_expon-a.expon())>51){ //limit of epsilon difference
		return *this;
	}
	_mantissa+=a.mantissa()*pow(2.0,(a.expon()-_expon)*1.0);
	fixParams();
	return *this;
}

inline doubleRepMantisa operator+(const doubleRepMantisa& a, const doubleRepMantisa& b){
	doubleRepMantisa temp(a);
	temp+=b;
	return temp;
}

inline doubleRepMantisa& doubleRepMantisa::operator-=(const doubleRepMantisa& a){
	doubleRepMantisa b(-a.mantissa(),a.expon());
	doubleRepMantisa me(_mantissa,_expon);
	me+=b;
	_mantissa=me.mantissa();
	_expon=me.expon();
	return *this;
}

inline doubleRepMantisa operator-(const doubleRepMantisa& a, const doubleRepMantisa& b){
	doubleRepMantisa temp(a);
	temp-=b;
	return temp;
}

inline doubleRepMantisa operator-(const doubleRepMantisa& a) {
        return doubleRepMantisa(0) - a;
}

inline doubleRepMantisa& doubleRepMantisa::operator*=(const doubleRepMantisa& a){
	_mantissa*=a.mantissa();
	_expon+=a.expon();
	fixParams();
	return *this;
}

inline doubleRepMantisa operator*(const doubleRepMantisa& a, const doubleRepMantisa& b){
	doubleRepMantisa temp(a);
	temp*=b;
	return temp;
}

inline doubleRepMantisa& doubleRepMantisa::operator/=(const doubleRepMantisa& a){
	_mantissa/=a.mantissa();
	_expon-=a.expon();
	fixParams();
	return *this;
}

inline doubleRepMantisa operator/(const doubleRepMantisa& a, const doubleRepMantisa& b){
	doubleRepMantisa temp(a);
	temp/=b;
	return temp;
}

/************************
 * Comparison operators *
 ************************/
inline bool operator==(const doubleRepMantisa& a, const doubleRepMantisa& b){
	return (a._mantissa==b._mantissa && a._expon==b._expon);
}
inline bool operator!=(const doubleRepMantisa& a, const doubleRepMantisa& b){
	return !(a==b);
}

inline bool operator<(const doubleRepMantisa& a, const doubleRepMantisa& b){
	// if the numbers have opposite signs
    if (a._mantissa*b._mantissa<0.0){
		if (a._mantissa<b._mantissa) {return true;}
		else {return false;}
    }
	// if the expon values are different
	if (a._expon!=b._expon) {
		// special case where one number is zero
		if (a._mantissa == 0.0) {
			if (b._mantissa > 0.0) {return true;}
			else {return false;}
		}
		if (b._mantissa == 0.0) {
			if (a._mantissa < 0.0) {return true;}
			else {return false;}
		}

		if (a._expon<b._expon) {
			if (a._mantissa > 0.0) {return true;}
			else {return false;}
		} else {
			if (a._mantissa < 0.0) {return true;}
			else {return false;}
		}
		// expon values are identical
	} else {
		return (a._mantissa < b._mantissa);
	}
}

inline bool operator>(const doubleRepMantisa& a, const doubleRepMantisa& b){
	// if the numbers have opposite signs
    if (a._mantissa*b._mantissa<0.0){
		if (a._mantissa>b._mantissa) {return true;}
		else {return false;}
    }
	// if the expon values are different
	if (a._expon!=b._expon) {
		// special case where one number is zero
		if (a._mantissa == 0.0) {
			if (b._mantissa < 0.0) {return true;}
			else {return false;}
		}
		if (b._mantissa == 0.0) {
			if (a._mantissa > 0.0) {return true;}
			else {return false;}
		}

		if (a._expon>b._expon) {
			if (a._mantissa > 0.0) {return true;}
			else {return false;}
		} else {
			if (a._mantissa < 0.0) {return true;}
			else {return false;}
		}
		// expon values are identical
	} else {
		return (a._mantissa > b._mantissa);
	}
}

inline bool operator<=(const doubleRepMantisa& a, const doubleRepMantisa& b){
	return !(a>b);
}

inline bool operator>=(const doubleRepMantisa& a, const doubleRepMantisa& b){
	return !(a<b);
}




ostream& operator<<(ostream &out, const doubleRepMantisa& a);
istream& operator>>(istream &in, doubleRepMantisa& a);

inline MDOUBLE log(const doubleRepMantisa& d) {return d.d_log();}

inline ostream &operator<<(ostream &out, const VdoubleRepMantisa &v){
  for (int j=0;j<v.size();++j)
    out<< v[j]<<" ";
  out <<endl;
  return(out);
}

inline ostream &operator<<(ostream &out, const VVdoubleRepMantisa &m){
  for (int i=0;i<m.size();++i)
    out<<m[i];
  out <<endl;
  return(out);
}

inline doubleRepMantisa pow(const doubleRepMantisa& d1, const doubleRepMantisa& d2) {
  return doubleRepMantisa(pow(convert(d1), convert(d2)));
}

inline doubleRepMantisa abs(const doubleRepMantisa& d) { 
  return doubleRepMantisa(abs(d._mantissa), d._expon);
}

inline doubleRepMantisa fabs(const doubleRepMantisa& d) { 
  return abs(d);
}

inline doubleRepMantisa exp(const doubleRepMantisa& d) {
  return doubleRepMantisa(exp(convert(d)));
}

inline doubleRepMantisa sqrt(const doubleRepMantisa& d) {
  return doubleRepMantisa(sqrt(convert(d)));
}





//inline const MDOUBLE convert (const MDOUBLE d) const  {return(d);}

#endif
#endif
