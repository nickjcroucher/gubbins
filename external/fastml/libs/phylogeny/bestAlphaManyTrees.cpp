// $Id: bestAlphaManyTrees.cpp 962 2006-11-07 15:13:34Z privmane $

// version 1.00
// last modified 3 Nov 2002

#include "bestAlphaManyTrees.h"
#include "bestAlpha.h"
#include "numRec.h"
#include "bblEMProportional.h"
#include "bblEMSeperate.h"
#include "logFile.h"
#include <iostream>
using namespace std;

#ifndef VERBOS
#define VERBOS
#endif


void bestAlpha::optimizeAlphaNG_EM_PROP(tree& et,
									vector<sequenceContainer>& sc,
									vector<stochasticProcess>& sp,
									const vector<Vdouble *> * weights,
									MDOUBLE & bestAlpha,
									MDOUBLE & likelihoodScore,
									const int maxIterations,
									const MDOUBLE epsilon){

	//LOG(5,<<" 1. bestAlpha::findBestAlpha"<<endl);
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE ax,bx,cx; // left, midle, right limit on alpha
	bx=1.5;	// the limits are becoming more narrow with time.
	ax=0;
	cx=5.0;
	MDOUBLE tol=0.01f;
	MDOUBLE bestA=0;
	int i;
	const int maxIterationsThisF = 50;
	for (i=0; i < maxIterationsThisF; ++i) {

		bblEMProportional bblEMprop1(et,sc,sp,weights,maxIterations,epsilon);
		MDOUBLE newL = bblEMprop1.getTreeLikelihood();
		
#ifdef VERBOS
		LOG(5,<<"Before optimizing alpha, L = "<<newL<<endl);
#endif

		MDOUBLE likeAfterAlphaOpt = -brent(ax,bx,cx, // NEW MINUS. CHECK
			C_evalAlphaManyTrees(et,sc,sp,weights),
			tol,
			&bestA); // THIS FUNCTION CHANGE SP, BUT YET ONE HAVE TO INSERT THE BEST ALPHAS.
		for (int z=0; z < sp.size();++z) {
			(static_cast<gammaDistribution*>(sp[z].distr()))->setAlpha(bestA);
		}

#ifdef VERBOS
		LOG(5,<<"After optimizing alpha, L = "<<likeAfterAlphaOpt<<endl);
		LOG(5,<<" best A = " << bestA<<endl);
#endif	
		newL = likeAfterAlphaOpt;
	
		

		if (newL > oldL+0.01) {
			oldL = newL;
		}
		else {
			if (newL > oldL) {
				likelihoodScore = newL;
				bestAlpha= bestA;
				return;
			}
			else {
				likelihoodScore = oldL;
				bestAlpha= bestA;
				return;
			}
		}
	}
	if (i == maxIterationsThisF) errorMsg::reportError(" to many iteration in function optimizeBranchLength");
}

/*
void findBestAlphaManyTrees::findBestAlphaFixedManyTrees(const vector<tree>& et,
					   vector<positionInfo>& pi,
					   const VVdouble * weights) {
	//LOG(5,<<" 1. bestAlpha::findBestAlpha"<<endl);
	MDOUBLE bestA=0;
	checkAllocation();
	MDOUBLE ax,bx,cx; // left, midle, right limit on alpha
	MDOUBLE tol;
	ax=0;bx=1.5;cx=2;
	tol=0.01f;
	_bestL = brent(ax,bx,cx,
			C_evalAlphaManyTrees(et,_pi,weights),
			tol,
			&bestA);
	_bestAlpha= bestA;
}

*/

void bestAlpha::optimizeAlphaNG_EM_SEP(
									vector<tree>& et,
									vector<sequenceContainer>& sc,
									vector<stochasticProcess>& sp,
									const vector<Vdouble *> * weights,
									MDOUBLE & bestAlpha,
									MDOUBLE & likelihoodScore,
									const int maxIterations,
									const MDOUBLE epsilon) {
	// SEPERATE ANALYSIS, 1 GAMMA
	//LOG(5,<<" 1. bestAlpha::findBestAlpha"<<endl);
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;
	MDOUBLE ax,bx,cx; // left, midle, right limit on alpha
	bx=1.5;	// the limits are becoming more narrow with time.
	ax=0;
	cx=5.0;
	MDOUBLE tol=0.01f;
	MDOUBLE bestA=0;
	const int maxIterationsThisF = 50;
	for (int i=0; i < maxIterationsThisF; ++i) {
		newL=0;
		LOG(3,<<"starting iteration "<<i<<endl);
		bblEMSeperate bblEMsep1(et,
									sc,
									sp,
									weights,
									maxIterations,
									epsilon);
		newL =bblEMsep1.getTreeLikelihood();			
#ifdef VERBOS
		LOG(5,<<"Before optimizing alpha, L = "<<newL<<endl);
#endif
		//MDOUBLE alphaB4optimizing = (static_cast<gammaDistribution*>(sp[0].distr()))->getAlpha();
		MDOUBLE likeAfterAlphaOpt = -brent(ax,bx,cx, // NEW MINUS - CHECK!
			C_evalAlphaManyTreesSep(et,sc,sp,weights),
			tol,
			&bestA);
		
		if (likeAfterAlphaOpt>newL) {
			for (int i=0; i < sc.size();++i) {
				(static_cast<gammaDistribution*>(sp[0].distr()))->setAlpha(bestA);
			}
			newL = likeAfterAlphaOpt;
		}
#ifdef VERBOS
		LOG(5,<<"After optimizing alpha, L = "<<newL<<endl);
#endif
		if (newL > oldL+0.01) {
			oldL = newL;
		}
		else {
			if (newL > oldL) {
				likelihoodScore = newL;
				bestAlpha= bestA;
				return;
			}
			else {
				likelihoodScore = oldL;
				bestAlpha= bestA;
				return;
			}
		}
	}
	errorMsg::reportError(" to many iteration in function optimizeBranchLength");
}

//==================== optimizing n alphas ==============================

void bestAlpha::optimizeAlphaNG_EM_PROP_n_alpha(tree& et,
									vector<sequenceContainer>& sc,
									vector<stochasticProcess>& sp,
									const vector<Vdouble *> * weights,
									vector<MDOUBLE> & bestAlphas,
									MDOUBLE & likelihoodScore,
									const int maxIterations,
									const MDOUBLE epsilon){

	//LOG(5,<<" 1. bestAlpha::findBestAlpha"<<endl);
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;
	MDOUBLE ax,bx,cx; // left, midle, right limit on alpha
	bx=1.5;	// the limits are becoming more narrow with time.
	ax=0;
	cx=5.0;
	vector<MDOUBLE> bestAs= bestAlphas;
	vector<MDOUBLE> newAlphas(sc.size(),0);
	int i;
	const int maxIterationsThisF = 50;
	for (i=0; i < maxIterationsThisF; ++i) {
#ifdef VERBOS
		LOG(5,<<" ============================ optimizing bbl (fixed alphas) ================= \n");
#endif
		newL=0;
		bblEMProportional bblem1(et,sc,sp,weights,maxIterations,epsilon);
		MDOUBLE tmpX =bblem1.getTreeLikelihood();
				
#ifdef VERBOS
		LOG(5,<<"likelihood of trees (sum)= "<<tmpX<<endl);
#endif
		newL =tmpX;
#ifdef VERBOS
		LOG(5,<<"Before optimizing alpha, L = "<<newL<<endl);
		LOG(5,<<" ============================ optimizing alphas ================= \n");
#endif
		const MDOUBLE upperBoundOnAlpha = 5;
		MDOUBLE likeAfterAlphaOpt = 0;
		for (int treeNumber =0; treeNumber<sc.size();++treeNumber) {
			bestAlphaFixedTree bestAlphaFixedTree1(et,
				sc[treeNumber],
				sp[treeNumber],
				weights?(*weights)[treeNumber]:NULL,
				upperBoundOnAlpha,
				epsilon);
			MDOUBLE tmpX = bestAlphaFixedTree1.getBestL();
#ifdef VERBOS
			LOG(5,<<"likelihood of tree "<<treeNumber<<" = "<<tmpX<<endl);
#endif
			newAlphas[treeNumber] = bestAlphaFixedTree1.getBestAlpha();
#ifdef VERBOS
			LOG(5,<<" best alpha tree number: "<<treeNumber<<" = "<<newAlphas[treeNumber]<<endl);
#endif
			likeAfterAlphaOpt +=tmpX;
		}
		 

		if (likeAfterAlphaOpt>newL) {
			for (int z=0; z < sp.size();++z) {
				(static_cast<gammaDistribution*>(sp[z].distr()))->setAlpha(newAlphas[z]);
			}
			newL = likeAfterAlphaOpt;
			bestAs = newAlphas;
		}
		
	#ifdef VERBOS
		LOG(5,<<"After optimizing alpha, L = "<<newL<<endl);
	#endif

		if (newL > oldL+0.01) {
			oldL = newL;
		}
		else {
			if (newL > oldL) {
				likelihoodScore = newL;
				bestAlphas= bestAs;
				return;
			}
			else {
				likelihoodScore = oldL;
				bestAlphas= bestAs;
				return;
			}
		}
	}
	if (i == maxIterationsThisF) {
		errorMsg::reportError(" to many iteration in function optimizeBranchLength");
	}
}

		//// CHECK:
		//MDOUBLE check_sum=0;
		//for (int k=0; k < sp.size(); ++k) {
		//	MDOUBLE check = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(et,sc[k],sp[k]);
		//	LOG(5,<<" CHECK = "<< check<<endl);
		//	check_sum+=check;
		//}
		//LOG(5,<<" check-sum = "<<check_sum<<endl);
		//// END CHECK
