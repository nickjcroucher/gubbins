#ifndef _MULTIPLE_STOCHASTIC_PROCESS
#define _MULTIPLE_STOCHASTIC_PROCESS

#include "stochasticProcess.h"


class multipleStochasticProcess {
public:
	multipleStochasticProcess();
	virtual ~multipleStochasticProcess();
	virtual MDOUBLE getProb(int spPlace) const;
	virtual stochasticProcess* getSp(int spPlace);
	virtual int getSPVecSize() const {return _spVec.size();}
	virtual void setSpVec(vector<stochasticProcess>& spVec);

    
protected:
	virtual void copy(const multipleStochasticProcess * pOther);
protected:
	vector<stochasticProcess> _spVec;
	Vdouble _spProb;
};
#endif
