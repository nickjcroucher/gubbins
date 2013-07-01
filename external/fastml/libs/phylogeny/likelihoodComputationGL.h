#ifndef ___LIKELIHOOD_COMPUTATION_GL
#define ___LIKELIHOOD_COMPUTATION_GL

#include "definitions.h"
#include "computePijComponent.h"
#include "sequenceContainer.h"
#include "suffStatComponent.h"
#include "unObservableData.h"
#include "computeUpAlg.h"


namespace likelihoodComputationGL {


	MDOUBLE getTreeLikelihoodAllPosAlphTheSame(const tree& tr,
		const sequenceContainer& sc,
		const vector<vector<stochasticProcess*> >& spVVec,
		const distribution * distGain, const distribution * distLoss,
		unObservableData *unObservableData_p);
	void fillPijAndUp(const tree& tr,
		const sequenceContainer& sc,
		const vector<vector<stochasticProcess*> >& spVVec,
		const distribution * distGain, const distribution * distLoss,
		vector<computePijGam>& pi_vec,
		vector<suffStatGlobalGam>& ssc_vec,
		vector<computeUpAlg>& cup_vec);
	MDOUBLE getProbOfPosUpIsFilledSelectionGam(const int pos,const tree& tr,
		const sequenceContainer& sc,
		const vector<vector<stochasticProcess*> >& spVVec, // only needed for sp.freq(let)
		const suffStatGlobalGamPos& cup,
		const distribution * distGain, const distribution * distLoss);
	
	
	MDOUBLE getTreeLikelihoodFromUp2(const tree& tr,
		const sequenceContainer& sc,
		const vector<vector<stochasticProcess*> >& spVVec,// only needed for sp.freq(let)
		const suffStatGlobalGam& cup,
		const distribution * distGain, const distribution * distLoss,unObservableData *unObservableData_p,
		Vdouble* posLike =NULL,
		const Vdouble * weights =NULL);
	MDOUBLE getTreeLikelihoodFromUp2(const tree& tr,
		const sequenceContainer& sc,
		const vector<vector<stochasticProcess*> >& spVVec,// only needed for sp.freq(let)
		const vector<suffStatGlobalGam>& cup_vec,
		const distribution * distGain, const distribution * distLoss,unObservableData *unObservableData_p,
		Vdouble* posLike =NULL,
		const Vdouble * weights =NULL);

// Error
	//MDOUBLE getTreeLikelihoodAllPosAlphTheSameNoComputeUp(const tree& tr,
	//	const sequenceContainer& sc,
	//	const vector<vector<stochasticProcess*> >& spVVec,
	//	const distribution * distGain, const distribution * distLoss,
	//	unObservableData *unObservableData_p);


///********************************************************************************************
//un-obervable data
//*********************************************************************************************/
//// used to fill the likelihood for the unobservable for each category
//	doubleRep getLofPos(const int pos,
//		const tree& tr,
//		const sequenceContainer& sc,
//		const computePijGam& pi,
//		const stochasticProcess& sp,
//		Vdouble& likePerCat);	// all the likdelhoodsPerCat and rateProb are filled
//// likelihood computation - full data (1)
//	MDOUBLE getTreeLikelihoodAllPosAlphTheSame(const tree& tr,
//		const sequenceContainer& sc,
//		const stochasticProcess& sp,
//		const Vdouble * const weights,
//		Vdouble *pLforMissingDataPerCat=NULL);
////	likelihood computation - per pos (1.1)
//	doubleRep getLofPos(const int pos,					// this function is used
//		const tree& tr,					// when gamma, and the br-len
//		const sequenceContainer& sc,		// are the same for all pos.
//		const computePijGam& pi,
//		const stochasticProcess& sp,
//		Vdouble *pLforMissingDataPerCat=NULL);
//// likelihood computation - per pos, per cat (1.1.1)
//	doubleRep getLofPos(const int pos,					// this function is used
//		const tree& tr,					// when the br-len
//		const sequenceContainer& sc,		// are the same for all
//		const computePijHom& pi,			// positions.
//		const stochasticProcess& sp);
//	
//	Vdouble getLofPosPerCat(const int pos,				// used when the likelihood given each category is needed, not only the sum
//		const tree& tr,
//		const sequenceContainer& sc,
//		const computePijGam& pi,
//		const stochasticProcess& sp);



};

#endif
