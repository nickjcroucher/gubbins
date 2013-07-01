// $Id: searchStatus.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___SEARCH_STATUS
#define ___SEARCH_STATUS

#include "definitions.h"

class searchStatus {
public:
	explicit searchStatus(const MDOUBLE startingTmp,const MDOUBLE factor);
	explicit searchStatus(){};
	void setParameters(const MDOUBLE tmp, const MDOUBLE factor) {
		_currentTmp=tmp;
		_factor=factor;
	}

	void tmpUp1(){_currentTmp *= _factor;}
	void tmpDown1(){_currentTmp /= _factor;}
	const MDOUBLE getTmp() const {return _currentTmp;}
	void setTmp(const MDOUBLE newTmp) {_currentTmp=newTmp;}
	virtual ~searchStatus(){}

private:
	MDOUBLE _currentTmp;
	MDOUBLE _factor;
};

#endif 


