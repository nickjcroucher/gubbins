// $Id: siteSpecificRate.h 4742 2008-08-19 17:40:56Z cohenofi $

#ifndef ___SITE_SPECIFIC_RATE
#define ___SITE_SPECIFIC_RATE

#include "definitions.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "computePijComponent.h"
#include "unObservableData.h"


// the function returns the total log-likelihood of the rates.
// it is used for computing the rates, when there is one tree common to 
// all positions and 1 stochastic process common to all position.

MDOUBLE computeML_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & likelihoodsV,
								   const sequenceContainer& sd,
								   const stochasticProcess& sp,
								   const tree& et,
								   const MDOUBLE maxRate=20.0f,
								   const MDOUBLE tol=0.0001f);

// this function is the same as the one above, but here, each site can have its
//own tree, or its own stochastic process.
//etVec: a vector of possible trees.
//spVec: a vector of possible stochastic processes.
//treeAttributesVec: defines which tree is assigned to a specific position.
//NOTE: the possible attributes are 1,2..., so that the tree for position i
//is etVec[treeAttributesVec[i]-1]
//The same is true for the stochastic process atributes vector.
MDOUBLE computeML_siteSpecificRate(Vdouble & ratesV,
						Vdouble & likelihoodsV,
						const Vint& spAttributesVec,
						const Vint& treeAttributesVec,
						const vector<tree> & etVec,
						const vector<const stochasticProcess *> & spVec,
						const sequenceContainer& sc,
						const MDOUBLE maxRate,
						const MDOUBLE tol);

// this function is the same as the one above, but here,
// there are only tree attributes.
MDOUBLE computeML_siteSpecificRate(Vdouble & ratesV,
						Vdouble & likelihoodsV,
						const Vint& treeAttributesVec,
						const vector<tree> & etVec,
						const stochasticProcess& sp,
						const sequenceContainer& sc,
						const MDOUBLE maxRate,
						const MDOUBLE tol);

// this function is the same as the one above, but here,
// there are only stochastic process attributes.
MDOUBLE computeML_siteSpecificRate(Vdouble & ratesV,
						Vdouble & likelihoodsV,
						const Vint& spAttributesVec,
						const tree & et,
						const vector<const stochasticProcess* > & spVec,
						const sequenceContainer& sc,
						const MDOUBLE maxRate,
						const MDOUBLE tol);

void computeML_siteSpecificRate(int pos,
								 const sequenceContainer& sc,
								 const stochasticProcess& sp,
								 const tree &et,
								 MDOUBLE& bestRate,
								 MDOUBLE& posL,
								 const MDOUBLE maxRate,
								 const MDOUBLE tol);

// BAYESIAN PART

// 1 sequence container, 1 tree, 1 position
void computeEB_EXP_siteSpecificRate(int pos,
									const sequenceContainer& sc,
									const stochasticProcess& sp,
									const computePijGam& cpg,
									const tree &et,
									MDOUBLE& bestRate,
									MDOUBLE & stdRate,
									MDOUBLE & lowerConf,
									MDOUBLE & upperConf,
									const MDOUBLE alphaConf,
									VVdouble* LpostPerCat=NULL,
									unObservableData* unObservableData_p=NULL);

// 1 stochastic process, 1 tree, all positions
void computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
									Vdouble & stdV,
									Vdouble & lowerBoundV,
									Vdouble & upperBoundV,
									const sequenceContainer& sc,
									const stochasticProcess& sp,
									const tree& et,
									const MDOUBLE alphaConf,
									VVdouble* LpostPerCat=NULL,
									unObservableData* unObservableData_p=NULL);


// many stochastic process, many tree, all positions
void computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & stdV,
								   Vdouble & lowerBoundV,
								   Vdouble & upperBoundV,
								   const Vint& spAttributesVec,
								   const Vint& treeAttributesVec,
							       const sequenceContainer& sc,
								   const vector<tree> & etVec,
								   const vector<const stochasticProcess *> & spVec,
								   const MDOUBLE alphaConf);

// many stochastic process, 1 tree, all positions
void computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & stdV,
								   Vdouble & lowerBoundV,
								   Vdouble & upperBoundV,
								   const Vint& spAttributesVec,
							       const sequenceContainer& sc,
								   const tree & et,
								   const vector<const stochasticProcess *> & spVec,
								   const MDOUBLE alphaConf);

// 1 stochastic process, many tree, all positions
void computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & stdV,
								   Vdouble & lowerBoundV,
								   Vdouble & upperBoundV,
								   const Vint& treeAttributesVec,
							       const sequenceContainer& sc,
								   const vector<tree> & etVec,
								   const stochasticProcess & sp,
								   const MDOUBLE alphaConf);
#endif

