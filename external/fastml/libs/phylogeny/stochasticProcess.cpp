// $Id: stochasticProcess.cpp 4660 2008-08-12 14:31:38Z cohenofi $

#include "stochasticProcess.h"
#include "errorMsg.h"

stochasticProcess& stochasticProcess::operator=(const stochasticProcess &otherStoc) {
	if (this != &otherStoc) {              // Check for self-assignment
		if (_pijAccelerator) delete _pijAccelerator;
		if (otherStoc._pijAccelerator)
		{
			pijAccelerator* p2 = otherStoc._pijAccelerator->clone();   // Create the new one FIRST...
			_pijAccelerator = p2;
		}
		else
			_pijAccelerator = NULL;

		if (_distr) delete _distr;
		if (otherStoc._distr)
		{
			distribution* d2 =  otherStoc._distr->clone();
			_distr = d2;
		}
		else{ 
			_distr = NULL;
			_isReversible = otherStoc.isReversible();
		}
	}
//	if (_distr) delete _distr;
//	_distr = new distribution(*otherStoc._distr);
    return *this;
}
   
	
stochasticProcess::stochasticProcess(const distribution *in_distr,const pijAccelerator *pijAccelerator, bool isReversible) :
	 _distr(in_distr->clone()), _pijAccelerator(pijAccelerator->clone()), _isReversible(isReversible){
	
}
	
stochasticProcess::stochasticProcess(const stochasticProcess& other):
	 _distr(NULL), _pijAccelerator(NULL){
		if (other._pijAccelerator != NULL) _pijAccelerator = other._pijAccelerator->clone();
		if (other._distr != NULL) _distr = other._distr->clone();
		_isReversible = other.isReversible();
}
	
stochasticProcess::~stochasticProcess() {
		delete _distr;
		delete _pijAccelerator;
}


void stochasticProcess::setDistribution(const distribution* in_distr)
{
	if (_distr)	delete _distr;
	if (in_distr == NULL) _distr = NULL;
	else _distr = in_distr->clone();
}
