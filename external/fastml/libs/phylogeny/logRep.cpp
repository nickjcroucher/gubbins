#ifdef LOGREP
#include "logRep.h"
#include <cmath>

//logRep::logRep()
//{
//	_log = VERYSMALL2;
//}

//logRep::logRep(MDOUBLE a){
//	_log = ((a==0.0) ? VERYSMALL2 : log(a));
//}


//logRep::logRep(const logRep& other): _log(other._log) {}



MDOUBLE convert(const logRep& a){
	return exp(a.getLog());
}




ostream& operator<<(ostream &out, const logRep& a){
	a.output(out);
    return out;
}
#endif
