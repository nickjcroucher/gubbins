#ifndef _NY_CODON_MODEL
#define _NY_CODON_MODEL

#include "replacementModel.h"
#include "fromQtoPt.h"
#include "codon.h"
#include "sequenceContainer.h"

class nyCodonModel : public replacementModel {
public:
	
	explicit nyCodonModel(const codon &inCodonAlpa, const Vdouble & codonFreq, bool bFillQ2pMatrix, const MDOUBLE synRate=1.0, const MDOUBLE nonsynRate=1.0, const MDOUBLE kappa=1.0);
	virtual ~nyCodonModel();
	const int alphabetSize() const {return inCodonAlpa.size();}
	virtual replacementModel* clone() const { return new nyCodonModel(*this); }
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
	
	MDOUBLE getKappa() const{return _kappa;}
	MDOUBLE getSynRate()const {return _synRate;}
	MDOUBLE getNonsynRate() const {return _nonsynRate;}
	VVdouble& getQ()  {return _Q;}

	void setKappa(const MDOUBLE k) {setParams(_synRate, _nonsynRate, k);}
	void setSynRate(const MDOUBLE synRate) {setParams(synRate, _nonsynRate, _kappa);}
	void setNonsynRate(const MDOUBLE nonsynRate) {setParams(_synRate, nonsynRate, _kappa);}
	void setParams(const MDOUBLE synRate, const MDOUBLE nonsynRate, const MDOUBLE kappa); 

	MDOUBLE getQij(const int i,const int j)const {return _Q[i][j];}
	void norm(MDOUBLE scale);
	MDOUBLE sumPijQij();
	bool isFillingQ2pMatrix() const {return _bFillQ2pMatrix;};
    void updateQ();

private:	
	void init(const codon &inCodonAlpa);
	void homogenousFreq(const codon &inCodonAlpa);
	//For each codon position calculates the frequency of each type of codon
	void calcNucFrequenciesPerPosition(const sequenceContainer& scNuc, const codon &inCodonAlpa);
	MDOUBLE getModelFreq(int fromCodon, int targetCodon); //this is not the codon frequency but the PI that is used in the Q matrix
	void initCodonFrequencies(const codon &co);

private:
	Vdouble _freq; //holds the fequncies of codons
	MDOUBLE _synRate; //syn
	MDOUBLE _nonsynRate; // nonsyn
	MDOUBLE _kappa; //tr/tv
	q2pt _q2pt;	
	VVdouble _Q;
	bool _bFillQ2pMatrix; // in case this model is part of a "father model" (multiple stochastic process) then 
						 //we don't want to fill the Q2P matrix every time we change a parameter.
						//The father model is responsible to normalize all its stochastic processes together and then fill the Q2P
	
};


#endif
