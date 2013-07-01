// $Id: NNiSep.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___NNI_SEP
#define ___NNI_SEP

#include "definitions.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "definitions.h"
#include "stochasticProcess.h"
#include <vector>
using namespace std;

class NNiSep {
public:
	explicit NNiSep(vector<sequenceContainer>& sc,
					 vector<stochasticProcess>& sp,
					const vector<Vdouble *> * weights,
					vector<char>* nodeNotToSwap);

	vector<tree> NNIstep(vector<tree> et);
	MDOUBLE bestScore(){ return _bestScore;} 
	void setOfstream(ostream* out);

private:
	vector<char>* _nodeNotToSwap;
	vector<tree> _bestTrees;
	MDOUBLE _bestScore;
	vector<sequenceContainer>& _sc;
	vector<stochasticProcess>& _sp;
	const vector<Vdouble *> * _weights;

	MDOUBLE evalTrees(vector<tree>& et);
	tree NNIswap1(tree et,tree::nodeP mynode);
	tree NNIswap2(tree et,tree::nodeP mynode);
	int _treeEvaluated;
	ostream* _out;

};
#endif
