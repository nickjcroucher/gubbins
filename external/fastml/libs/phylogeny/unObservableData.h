#ifndef ___unObservableData___GL
#define ___unObservableData___GL

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "gainLossAlphabet.h"
#include "computePijComponent.h"

/********************************************************************************************
unObservableData
*********************************************************************************************/
class unObservableData{
public:
	explicit unObservableData(const sequenceContainer& sc,const stochasticProcess* sp ,const gainLossAlphabet alph, const int minNumOfOnes);
	unObservableData(const unObservableData& other); //const
	virtual ~unObservableData(){};
	virtual unObservableData* clone() const {return new unObservableData(*this);}
	Vdouble* getpLforMissingDataPerCat();
	Vdouble getLforMissingDataPerCat();
	MDOUBLE getlogLforMissingData();
	int getNumOfUnObservablePatterns();
	void setLforMissingData(const tree& _tr, const stochasticProcess* _sp);
	//void setLforMissingData(const tree& _tr, const stochasticProcess* _sp);
	void setLforMissingData(const tree& _tr, const vector<vector<stochasticProcess*> >& spVVec,	const distribution * distGain, const distribution* distLoss);



	//MDOUBLE getCorrectedLikelihood(MDOUBLE likePre){return }


protected:
//func

protected:
//members
	sequenceContainer _scZero;
	Vdouble _LforMissingDataPerCat;			// used foreach rate category
	MDOUBLE _logLforMissingData;
	computePijGam _pi;
};


#endif
