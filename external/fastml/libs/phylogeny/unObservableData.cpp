#include "unObservableData.h"
#include "likelihoodComputation.h"
#include "likelihoodComputationGL.h"
#include <math.h>


using namespace std;

unObservableData::unObservableData(const sequenceContainer& sc,const stochasticProcess* sp ,const gainLossAlphabet alph, const int minNumOfOnes)
{
	_scZero.startZeroSequenceContainerGL(sc,alph, minNumOfOnes);
	_LforMissingDataPerCat.resize(sp->categories());
}

unObservableData::unObservableData(const unObservableData& other) //const
{
	_scZero = other._scZero;
	_pi = other._pi;
	_logLforMissingData = other._logLforMissingData;
	_LforMissingDataPerCat = other._LforMissingDataPerCat;
}
Vdouble* unObservableData::getpLforMissingDataPerCat(){return &_LforMissingDataPerCat;}
Vdouble unObservableData::getLforMissingDataPerCat(){return _LforMissingDataPerCat;}
MDOUBLE unObservableData::getlogLforMissingData(){return _logLforMissingData;}
int unObservableData::getNumOfUnObservablePatterns(){return _scZero.seqLen();}


//void unObservableData::setLforMissingData(const tree& _tr, const stochasticProcess* _sp){
//	_pi.fillPij(_tr,*_sp);
//// NOTE: The "perCat" is out
//	_LforMissingDataPerCat = likelihoodComputation::getLofPosPerCat(0,_tr,_scZero,_pi,*_sp);	// L * sp.ratesProb(i)
//	_logLforMissingData = 0;
//	for (int i=0; i < _sp->categories();++i) {
//		_logLforMissingData += _LforMissingDataPerCat[i];
//	}
//	_logLforMissingData = log(_logLforMissingData);
//}

/********************************************************************************************
*********************************************************************************************/
void unObservableData::setLforMissingData(const tree& tr, const stochasticProcess* sp){
	_pi.fillPij(tr,*sp);
	_logLforMissingData = 0;
	for(int pos=0; pos<_scZero.seqLen(); ++pos){
		_logLforMissingData += convert(likelihoodComputation::getLofPos(pos,tr,_scZero,_pi,*sp));	
	}
	_logLforMissingData = log(_logLforMissingData);
}
/********************************************************************************************
*********************************************************************************************/
void unObservableData::setLforMissingData(const tree& tr, const vector<vector<stochasticProcess*> >& spVVec,	
						const distribution* distGain, const distribution* distLoss)
{
	
	_logLforMissingData = 0;
	int numOfRateCategories = spVVec[0][0]->categories();
	vector<computePijGam> pi_vec(numOfRateCategories);
	vector<suffStatGlobalGam> ssc_vec(numOfRateCategories);
	vector<computeUpAlg> cup_vec(numOfRateCategories);
	likelihoodComputationGL::fillPijAndUp(tr,_scZero, spVVec,distGain,distLoss,pi_vec,ssc_vec,cup_vec);
	
	for (int k=0; k < _scZero.seqLen(); ++k) {
		MDOUBLE resGivenRate = 0.0;
		MDOUBLE lnL = 0;
		for(int rateIndex=0 ; rateIndex<numOfRateCategories; ++rateIndex){
			lnL = log(likelihoodComputationGL::getProbOfPosUpIsFilledSelectionGam(k,//pos,
				tr,//const tree& 
				_scZero,// sequenceContainer& sc,
				spVVec,	// only needed for sp.freq(let)
				ssc_vec[rateIndex][k],//const computePijGam& ,
				distGain, distLoss)); // distributions
			resGivenRate += lnL * spVVec[0][0]->ratesProb(rateIndex);
		}
		_logLforMissingData += exp(resGivenRate);
	}
	_logLforMissingData = log(_logLforMissingData);
	//for(int rateIndex=0 ; rateIndex<numOfRateCategories; ++rateIndex){
	//	_logLforMissingData += likelihoodComputationGL::getTreeLikelihoodFromUp2(tr,_scZero,spVVec,ssc_vec[rateIndex], distGain,distLoss,NULL)
	//		* spVVec[0][0]->ratesProb(rateIndex);
	//}
}

