#ifndef __LOG_REP_H
#define __LOG_REP_H

#ifdef LOGREP

#include "definitions.h"
#include "AddLog.h"



#include <iostream>
#include <cmath>
using namespace std;

/* logRep: enables working with much larger or smaller numbers than normally possible
by the regular double representation
 * Representation of a number x by the log of x
 Note: Base is 2!!
 WARNING: Note that logRep can only be used for positive values 
 (such as probablities) - you can't have the log of a negative!
 For a general real number use class doubleRep.
 */ 

class logRep{
public:
	
	logRep() : _log(VERYSMALL){}
	logRep(MDOUBLE a) {_log = ((a==0.0) ? VERYSMALL : log(a));}
	logRep(const logRep& other) : _log(other._log) {}
	logRep* clone() {return new logRep(*this);}

	void output(ostream &out) const{ out<<exp(_log);}
    
	friend MDOUBLE convert(const logRep& a);
	//inline MDOUBLE convert();
	inline logRep& operator=(const logRep& a);
	inline logRep& operator+=(logRep a);
	friend inline logRep operator+(const logRep& a, const logRep& b);
	inline logRep& operator-=(const logRep& a);
	friend inline logRep operator-(const logRep& a, const logRep& b);
	inline logRep& operator*=(const logRep& a);
	friend inline logRep operator*(const logRep& a, const logRep& b);
	inline logRep& operator/=(const logRep& a);
	friend inline logRep operator/(const logRep& a, const logRep& b);

	friend inline bool operator==(const logRep& a, const logRep& b);
	friend inline bool operator!=(const logRep& a, const logRep& b);
	friend inline bool operator<(const logRep& a, const logRep& b);
	friend inline bool operator<=(const logRep& a, const logRep& b);
	friend inline bool operator>(const logRep& a, const logRep& b);
	friend inline bool operator>=(const logRep& a, const logRep& b);
	friend inline MDOUBLE log(const logRep& d);
	
private:
	const MDOUBLE getLog() const {return _log;}

private:
	MDOUBLE _log;
	//static tAddLog_Precompute _add;
	
};
	
inline logRep& logRep::operator=(const logRep& a){
	_log=a.getLog();
	return *this;
}

//inline MDOUBLE convert(){
//	return exp(_log);
//}

// Original version by Adi Stern
inline logRep& logRep::operator+=(logRep a){
	if (_log == VERYSMALL) 
		_log = a._log;
	else if (a._log == VERYSMALL ) return *this;
	else _log = AddLog(_log, a._log); 
	return *this;
}

inline logRep operator+(const logRep& a, const logRep& b){
	logRep temp(a);
	temp+=b;
	return temp;
}

inline logRep& logRep::operator*=(const logRep& a){
	if ((_log == VERYSMALL) || (a._log== VERYSMALL )){
		_log = VERYSMALL;
		return *this;
	}
	_log+=a._log;
	return *this;
}

inline logRep operator*(const logRep& a, const logRep& b){
	logRep temp(a);
	temp*=b;
	return temp;
}

inline logRep& logRep::operator/=(const logRep& a){
	_log-=a._log;
	return *this;
}

inline logRep operator/(const logRep& a, const logRep& b){
	logRep temp(a);
	temp/=b;
	return temp;
}

/************************
 * Comparison operators *
 ************************/
inline bool operator==(const logRep& a, const logRep& b){
	return (a.getLog()==b.getLog());
}
inline bool operator!=(const logRep& a, const logRep& b){
	return !(a==b);
}

inline bool operator<(const logRep& a, const logRep& b){
	if (a.getLog()<b.getLog()) {return true;}
	else {return false;}
    
}

inline bool operator>(const logRep& a, const logRep& b){
	
	if (a.getLog()>b.getLog()) {return true;}
	else {return false;}
   
}

inline bool operator<=(const logRep& a, const logRep& b){
	return !(a>b);
}

inline bool operator>=(const logRep& a, const logRep& b){
	return !(a<b);
}

ostream& operator<<(ostream &out, const logRep& a);

inline MDOUBLE log(const logRep& d) {return d.getLog();}

inline ostream &operator<<(ostream &out, const VlogRep &v){
  for (int j=0;j<v.size();++j)
    out<< v[j]<<" ";
  out <<endl;
  return(out);
}

inline ostream &operator<<(ostream &out, const VVlogRep &m){
  for (int i=0;i<m.size();++i)
    out<<m[i];
  out <<endl;
  return(out);
}
#endif
#endif
