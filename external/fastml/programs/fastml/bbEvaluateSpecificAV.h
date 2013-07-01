#if !defined ___BB__EVALUATE_SPECIFIC_AV__
#define ___BB__EVALUATE_SPECIFIC_AV__

#include "bb_options.h"
#include "computePijComponent.h" 
#include "suffStatComponent.h"
#include "sequence.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "tree.h"
#include "seqContainerTreeMap.h"

#include <vector>
using namespace std;

class bbEvaluateSpecificAV {

public:
	explicit	bbEvaluateSpecificAV(
					const tree& et,
					const stochasticProcess& sp,
					const sequenceContainer& sc,
					const computePijGam& cpij);
	virtual ~bbEvaluateSpecificAV();

	MDOUBLE evaluateSpecificAv(	const int pos,
								const vector<sequence>* ancestralSequences);
private:
	const tree& _et;
	const stochasticProcess& _sp;
	const computePijGam& _bbcpij;
	int _alphabetSize;
	int _pos;
	const sequenceContainer& _sc;
	seqContainerTreeMap * _sctm;


	const vector<sequence>* _ancss;

	MDOUBLE recursiveEvaluateSpecificAv(
					const int pos,
					const tree::nodeP thisNode);

	MDOUBLE recursiveEvaluateSpecificAv(const int pos,
											const tree::nodeP thisNode,
											const int categor);
	VVdouble _Lvec; // inodes * letter
	
};

#endif
