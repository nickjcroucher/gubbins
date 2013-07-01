#if !defined ___BB__FIND_BEST_AV_DYN_PROG
#define ___BB__FIND_BEST_AV_DYN_PROG


#include "bb_options.h"
#include "computePijComponent.h" 
#include "suffStatComponent.h"
#include "sequence.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "seqContainerTreeMap.h"

class bbfindBestAVDynProg {
public:
	explicit bbfindBestAVDynProg(const tree* et,
					const stochasticProcess *sp,
					const sequenceContainer& sc,
					const computePijGam* cpij);
	virtual ~bbfindBestAVDynProg();

	MDOUBLE evaluateSpecificAvDP(	const int pos,
									const vector<sequence>* ancestralSequences,
									const int rateCategory
	);

private:
	const tree* _et;
	const stochasticProcess* _sp;
	const computePijGam* _bbcpij;
	int _alphabetSize;
	int _pos;
	seqContainerTreeMap * _sctm;
	const sequenceContainer& _sc;

	const vector<sequence>* _ancss;

	void recursiveComputeLandC(	const int pos,
								const tree::nodeP inNode,
								const int rateCategor);
	VVdouble _jointLval; // inodes * letter
	VVdouble _jointCval; // inodes * letter
};

#endif
