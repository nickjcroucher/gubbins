// $Id: treeInference.cpp 962 2006-11-07 15:13:34Z privmane $
#include "treeInference.h"
#include "likeDist.h"
#include "distanceTable.h"

tree treeInference::computeNJtreeWithLikeDist(const stochasticProcess &sp, const sequenceContainer &sc, 
				   const tree * const constraintTreePtr, const vector<MDOUBLE> * const weights) {

	likeDist ld( sp, 0.01);
	VVdouble disTab;
	vector<string> vNames;
	giveDistanceTable(&ld,sc,disTab,vNames,weights);
	NJalg nj1;
	return (nj1.computeTree(disTab,vNames,constraintTreePtr));
}

