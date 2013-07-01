#ifndef ___SUFF_STAT_GAMMA_MIXTURE
#define ___SUFF_STAT_GAMMA_MIXTURE
/************************************************************
The suffStatGammaMixture class is used to obtain the sufficient statistics 
that are neccessary for the EM algorithm to compute the mixture distribution parameters.
The following notations are used below:
P(h[i]=k): the probability that position i belongs to the Kth Gamma component.
teta_t: the current mixture distribution parameters (the alpha, beta, and the probability of each component).

There are 3 sufficient statistics:
M_k:	the expected number of positions belong to the Kth component.
		sigma(i = 1 to seqLen){P(h[i] = k|data, cur_mixtureDistribution)} 
A_k:	sigma(i = 1 to seqLen){P(h[i] = k|data, cur_mixtureDistribution) * E[r|h[i] = k, data, cur_mixtureDistribution]}
B_k:	sigma(i = 1 to seqLen){P(h[i] = k|data, cur_mixtureDistribution) * E[log(r)|h[i] = k, data, cur_mixtureDistribution]}
************************************************************/
#include "definitions.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "tree.h"
#include "mixtureDistribution.h"
#include "computePijComponent.h"

class suffStatGammaMixture{

public:

	explicit suffStatGammaMixture(const stochasticProcess& cur_sp, const sequenceContainer& sc, const tree& inTree);
	virtual ~suffStatGammaMixture();

	void computeStatistics();

	void plotStatistics(ofstream & outF);
	MDOUBLE getMk(int comp) const {return _MkVec[comp];}
	MDOUBLE getAk(int comp) const {return _AkVec[comp];}
	MDOUBLE getBk(int comp) const {return _BkVec[comp];}
	MDOUBLE computeQ(); 
	MDOUBLE computeQ2(); 


private:
	MDOUBLE computeStatisticsForComponent(int pos, int componentNum, const computePijGam& cpg);
	void allocatePlaceForSuffStat();
	void computePijForEachComponent(vector<computePijGam>& cpgVec,vector<stochasticProcess>& spVec);

private:	
	Vdouble _MkVec;
	Vdouble _AkVec;
	Vdouble _BkVec;

	const stochasticProcess* _pSp;
	const sequenceContainer* _pSc;
	const tree* _pTree;
};



#endif

