#ifndef _W_YANG_MODEL
#define _W_YANG_MODEL

#include "replacementModel.h"
#include "fromQtoPt.h"
#include "codon.h"


class wYangModel : public replacementModel {
public:
	explicit wYangModel(const MDOUBLE inW, const MDOUBLE inK,bool globalW, codon * coAlpha);
	explicit wYangModel(const MDOUBLE inW, const MDOUBLE inK, const Vdouble& freq,bool globalW, codon *coAlpha);
	explicit wYangModel(const wYangModel &other): _coAlpha(NULL) {(*this) = other;}
	virtual wYangModel& operator=(const wYangModel &other);
	virtual wYangModel* clone() const { return new wYangModel(*this); }
	virtual ~wYangModel() {
		if (_coAlpha) 
			delete _coAlpha;
	}

	const int alphabetSize() const {return _freq.size();}
	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const {
		return _q2pt.Pij_t(i,j,d);
	}
	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const{
		return _q2pt.dPij_dt(i,j,d);
	}
	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const{
		return _q2pt.d2Pij_dt2(i,j,d);
	}
	const MDOUBLE freq(const int i) const {return _freq[i];};
	void setK(const MDOUBLE newK) { _k = newK; updateQ();}
	void setW(const MDOUBLE newW) { _w = newW;updateQ();}
	void homogenousFreq(){ _freq.erase(_freq.begin(),_freq.end()),_freq.resize(alphabetSize(),1.0/alphabetSize());}

	MDOUBLE getK() const {return _k;}
	MDOUBLE getW() const {return _w;}
	
	MDOUBLE getQij(const int i,const int j)const {return _Q[i][j];}
	void setGlobalW(bool globalW){_globalW = globalW;}
	void norm(MDOUBLE scale);
	MDOUBLE sumPijQij();
private:
	void updateQ();


private:
	
	MDOUBLE _w; //selection factor.
	MDOUBLE _k; // Tr/Tv ratio.
	q2pt _q2pt;	
	VVdouble _Q;
	bool _globalW;   //false when compute w per site
	Vdouble _freq;
	codon *_coAlpha;
};


#endif
