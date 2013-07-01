#ifndef _GTR_MODEL
#define _GTR_MODEL

#include "replacementModel.h"
#include "fromQtoPt.h"

class gtrModel : public replacementModel {
public:
	enum modelElements {a = 0,c,g,t};
	explicit gtrModel(const Vdouble& freq,
					  const MDOUBLE a2c = 0.25,
					  const MDOUBLE a2g = 0.25,
					  const MDOUBLE a2t = 0.25,
					  const MDOUBLE c2g = 0.25,
					  const MDOUBLE c2t = 0.25,
					  const MDOUBLE g2t = 0.25);
	virtual replacementModel* clone() const { return new gtrModel(*this); }
	virtual gtrModel& operator=(const gtrModel &other);
	explicit gtrModel(const gtrModel &other);
	const int alphabetSize() const {return _freq.size();}
	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const {return _q2pt.Pij_t(i,j,d);}
	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const{return _q2pt.dPij_dt(i,j,d);}
	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const{return _q2pt.d2Pij_dt2(i,j,d);}
	const MDOUBLE freq(const int i) const {return _freq[i];};
	void set_a2c(const MDOUBLE a2c);
	void set_a2g(const MDOUBLE a2g);
	void set_a2t(const MDOUBLE a2t);
	void set_c2g(const MDOUBLE c2g);
	void set_c2t(const MDOUBLE c2t);
	void set_g2t(const MDOUBLE g2t);
	MDOUBLE get_a2c() const;
	MDOUBLE get_a2g() const;
	MDOUBLE get_a2t() const;
	MDOUBLE get_c2g() const;
	MDOUBLE get_c2t() const;
	MDOUBLE get_g2t() const;
	const VVdouble& getQ() const {return _Q;}
    

private:
	void updateQ(const MDOUBLE a2c,const MDOUBLE a2g,const MDOUBLE a2t,const MDOUBLE c2g,const MDOUBLE c2t,const MDOUBLE g2t);
	void norm(const MDOUBLE scale);
	MDOUBLE sumPijQij();

private:
	VVdouble _Q;
	Vdouble _freq;
	q2pt _q2pt;	
	MDOUBLE _a2c;
	MDOUBLE _a2g;
	MDOUBLE _a2t;
	MDOUBLE _c2g;
	MDOUBLE _c2t;
	MDOUBLE _g2t;
};
#endif






