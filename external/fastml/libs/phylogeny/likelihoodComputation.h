// $Id: likelihoodComputation.h 5058 2008-10-19 15:55:24Z cohenofi $

#ifndef ___LIKELIHOOD_COMPUTATION
#define ___LIKELIHOOD_COMPUTATION

#include "definitions.h"
#include "computePijComponent.h"
#include "sequenceContainer.h"
#include "suffStatComponent.h"
#include "unObservableData.h"


namespace likelihoodComputation {
// likelihood computation - full data (1)
	MDOUBLE getTreeLikelihoodAllPosAlphTheSame(const tree& et,
		const sequenceContainer& sc,
		const stochasticProcess& sp,
		const Vdouble * const weights = NULL,
		unObservableData *unObservableData_p=NULL);
//	likelihood computation - per pos (1.1)
	doubleRep getLofPos(const int pos,					// this function is used
		const tree& et,					// when gamma, and the br-len
		const sequenceContainer& sc,		// are the same for all pos.
		const computePijGam& pi,
		const stochasticProcess& sp,
		unObservableData *unObservableData_p=NULL);
// likelihood computation - per pos, per cat (1.1.1)
	doubleRep getLofPos(const int pos,					// this function is used
		const tree& et,					// when the br-len
		const sequenceContainer& sc,		// are the same for all
		const computePijHom& pi,			// positions.
		const stochasticProcess& sp,
		unObservableData *unObservableData_p=NULL);



	// used when the likelihood given each category is needed, not only the sum	
	Vdouble getLofPosPerCat(const int pos,				
		const tree& et,
		const sequenceContainer& sc,
		const computePijGam& pi,
		const stochasticProcess& sp);
	// used to fill the likelihood for the unobservable for each category
	doubleRep getLofPos(const int pos,
		const tree& et,
		const sequenceContainer& sc,
		const computePijGam& pi,
		const stochasticProcess& sp,
		Vdouble& likePerCat);	// all the likdelhoodsPerCat and rateProb are filled




// --------------------------------------------------------------------------------
// this function should be used only when the branch lengths are not the same for
// all positions. Otherwise, computePijHom should be calculated once,
// and be used for all calls. In this function, computePijHom is being computed for
// each position.
doubleRep getLofPosHomModelEachSiteDifferentRate(const int pos,			
					  const tree& et,					
					  const sequenceContainer& sc,		
					  const stochasticProcess& sp);
// ---------------------------------------------------------------------------------


// --------------------------------------------------------------------------------
// this function should be used only when the branch lengths are not the same for
// all positions. Otherwise, computePijHom should be calculated once,
// and be used for all calls. In this function, computePijHom is being computed for
// each position.
doubleRep getLofPosGamModelEachSiteDifferentRate(const int pos,
					  const tree& et,
					  const sequenceContainer& sc,
					  const stochasticProcess& sp);
// --------------------------------------------------------------------------------


	doubleRep getLofPos(const int pos,					// with a site specific rate.
					  const tree& et,
					  const sequenceContainer& sc,
					  const stochasticProcess& sp,
					  const MDOUBLE gRate);
	doubleRep getProbOfPosWhenUpIsFilledHom(const int pos,	// to be used for homogenous model
					  const tree& et,
					  const sequenceContainer& sc,
					  const stochasticProcess& sp,
					  const suffStatGlobalHomPos& ssc);
	doubleRep getProbOfPosWhenUpIsFilledGam(const int pos, // to be used for Gamma model.
						const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGamPos& cup);

  doubleRep getLofPosAndPosteriorOfRates(const int pos,
						const tree& et,
						const sequenceContainer& sc,
						const computePijGam& pi,
						const stochasticProcess& sp,
						VdoubleRep& postrior);

	MDOUBLE getTreeLikelihoodFromUp(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGam& cup,
						const Vdouble * weights =0 );
	
	MDOUBLE getTreeLikelihoodFromUp2(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGam& cup,
						VdoubleRep& posLike, // fill this vector with each position likelihood but without the weights.
						const Vdouble * weights=0,
						unObservableData* unObservableData_p=NULL);
	// fill this vector with each position posterior rate (p(r|D))
	// but without the weights.
	// the weights are used only because we return the likelihood 
	// (this takes these weights into account).
	MDOUBLE getPosteriorOfRates(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGam& cup,
						VVdoubleRep& posteriorLike, 
						const Vdouble * weights = NULL);

	MDOUBLE getPosteriorOfRates(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						VVdoubleRep& posteriorLike, 
						const Vdouble * weights = NULL);

  // fill the posteriorLike matrix with each position posterior rate (p(r|D))
  // and the LLPP,  but without the weights.
  MDOUBLE getPosteriorOfRatesAndLLPP(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGam& cup,
						VVdoubleRep& posteriorLike, 
						VdoubleRep& LLPerPos, 
						const Vdouble * weights=NULL);
	// From Itay M.
	// this function forces non gamma computation of likelihoods from up.
	// i.e., even if the stochastic process is really gamma - the likelihood is computed as if there's no gamma.
	MDOUBLE getTreeLikelihoodFromUpSpecifcRates(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalHom& cup,
						VdoubleRep& posLike, // fill this vector with each position likelihood but without the weights.
						const Vdouble * weights = NULL);

	// added from main semphy on 23.5.2005 (eyal privman + matan ninio).
	MDOUBLE computeLikelihoodAndLikelihoodPerPosition(const sequenceContainer &sc, const tree &et, 
												  const stochasticProcess &sp, Vdouble &LLPerPos);
    MDOUBLE getTreeLikelihoodFromPosteriorAndAlpha(const MDOUBLE alpha,
						const Vdouble originalBounderi,
						const VVdouble& posteriorLike,
						const VdoubleRep& LLPP,
						const Vdouble* weights);
    


};



#endif

