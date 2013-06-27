// $Id: fastStartTree.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___FAST_START_TREE
#define ___FAST_START_TREE

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include <iostream>

using namespace std;



tree getBestMLTreeFromManyNJtrees(sequenceContainer & allTogether,
								stochasticProcess& sp,
								const int numOfNJtrees,
								const MDOUBLE tmpForStartingTreeSearch,
								const MDOUBLE epslionWeights,
								ostream& out);


#endif
