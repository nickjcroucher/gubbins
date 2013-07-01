#include "multipleStochasticProcess.h"
#include "errorMsg.h"

multipleStochasticProcess::multipleStochasticProcess()
{
}


multipleStochasticProcess::~multipleStochasticProcess()
{
}


void multipleStochasticProcess::copy(const multipleStochasticProcess *pOther)
{
	_spVec = pOther->_spVec;
	_spProb = pOther->_spProb;
}


MDOUBLE multipleStochasticProcess::getProb(int spPlace) const {
	if (spPlace >= _spProb.size())
		errorMsg::reportError("error in multipleStochasticProcess::getProb");
	return _spProb[spPlace];
}

stochasticProcess* multipleStochasticProcess::getSp(int spPlace) {
	if (spPlace >= _spVec.size())
		errorMsg::reportError("error in multipleStochasticProcess::getSp");
	return &_spVec[spPlace];
}

void multipleStochasticProcess::setSpVec(vector<stochasticProcess>& spVec)
{
	_spVec.clear();
	_spVec = spVec;
}

