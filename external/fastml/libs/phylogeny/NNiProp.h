// $Id: NNiProp.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___NNI_PROP
#define ___NNI_PROP
#include "definitions.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "definitions.h"
#include "stochasticProcess.h"
#include <vector>
using namespace std;

class NNiProp {
public:
	explicit NNiProp(vector<sequenceContainer>& sc,
					 vector<stochasticProcess>& sp,
					const vector<Vdouble *> * weights,
					vector<char>* nodeNotToSwap);

	tree NNIstep(tree et);
	MDOUBLE bestScore(){ return _bestScore;} 
	void setOfstream(ostream* out);
private:
	ostream* _out;
	vector<char> * _nodeNotToSwap;
private:
	tree _bestTree;
	MDOUBLE _bestScore;
	vector<sequenceContainer>& _sc;
	vector<stochasticProcess>& _sp;
	const vector<Vdouble *> * _weights;

	MDOUBLE evalTree(tree& et);
	tree NNIswap1(tree et,tree::nodeP mynode);
	tree NNIswap2(tree et,tree::nodeP mynode);
	int _treeEvaluated;
	
};
#endif
