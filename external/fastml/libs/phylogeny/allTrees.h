// $Id: allTrees.h 1731 2007-02-26 13:45:23Z itaymay $

#ifndef ___ALL_TREES
#define ___ALL_TREES

#include "definitions.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include <vector>
using namespace std;

void get3seqTreeAndIdLeftVec(const sequenceContainer* sc,
							 tree& starT,
							 vector<int>& idList);

tree getAnewTreeFrom(	const tree& et,
							tree::nodeP & mynode,
							vector<int> & idLeft,
							const string& nameToAdd);
class allTrees {
public:
	explicit allTrees(bool keepAllTrees = false);
	MDOUBLE getBestScore() {return _bestScore;}
	tree getBestTree() {return _bestTree;}

	void getAllTreesAndLikelihoods(vector<tree>& resTree,VdoubleRep & scores) {
		resTree = _allPossibleTrees;
		scores = _allPossibleScores;
	}

	void recursiveFind(	tree et,
						const stochasticProcess& sp,
						const sequenceContainer& sc,
						vector<int> idLeft,
						const Vdouble * weights = NULL,
						const int maxIterations=1000,
						const MDOUBLE epsilon=0.05);

	void recursiveFind(	const sequenceContainer* sc,
						const stochasticProcess* sp,
						const Vdouble * weights = NULL,
						const int maxIterations=1000,
						const MDOUBLE epsilon=0.05); // one tree.
	


private:
	tree _bestTree;
	MDOUBLE _bestScore;
	vector<tree> _allPossibleTrees;
	vector<doubleRep> _allPossibleScores;
	const bool _keepAllTrees;


	MDOUBLE evalTree(tree& et,
					const stochasticProcess& sp,
					const sequenceContainer& sc,
					const int maxIterations,
					const MDOUBLE epsilon,
					const Vdouble * weights = NULL);




};
#endif

