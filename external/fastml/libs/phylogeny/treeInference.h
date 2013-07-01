// $Id: treeInference.h 962 2006-11-07 15:13:34Z privmane $
// 

// version 1.01
// last modified 23 May 2005

#ifndef ___TREE_INFERENCE
#define ___TREE_INFERENCE

#include "definitions.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "nj.h"
#include <vector>
using namespace std;

class treeInference {
public:
	static tree computeNJtreeWithLikeDist(const stochasticProcess &sp, const sequenceContainer &sc, 
				   const tree * const constraintTreePtr = NULL, const vector<MDOUBLE> * const weights = NULL);

};
#endif


