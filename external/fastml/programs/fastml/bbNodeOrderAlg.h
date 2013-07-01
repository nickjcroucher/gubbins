#if !defined ___BB__NODE_ORDER_ALG__
#define ___BB__NODE_ORDER_ALG__

#include "definitions.h"
#include "bb_options.h"
#include "computePijComponent.h" 
#include "suffStatComponent.h"
#include "sequence.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"

class bbNodeOrderAlg {
public:
	explicit bbNodeOrderAlg(const tree& et,
					const stochasticProcess &sp,
					const sequenceContainer& sc,
					const computePijGam& cpij,
					const MDOUBLE computeAgainExactTreshold);
	virtual ~bbNodeOrderAlg();
	void getNextNodeAndCharOrder(tree::nodeP &nextNode,
								 vector<int> &charOrder,
								 vector<sequence> &ancestralSequences,
								 const int pos,
								 const bool firstTime,
								 doubleRep& exactVal);
	void putBack(tree::nodeP& node2check,const doubleRep & exactVal);

private:
	const tree& _et;
	const stochasticProcess& _sp;
	const computePijGam& _cpij;
	const sequenceContainer& _sc;
	suffStatGlobalGamPos cmarginalbb;
	suffStatGlobalGamPos cupbb;
	suffStatGlobalGamPos cdownbb;

	MDOUBLE _computeAgainExactTreshold;
	int _alphabetSize;
	int _pos;
	vector<tree::nodeP> _nodesLeft;
	vector<doubleRep> _nodesLeftExact;

	void recalculateExact(	vector<sequence> &ancestralSequences,
							const int pos);
	vector<int> findBestOrderInNode(const tree::nodeP node2check,
												const int pos) const;
	void rankRemainingNodesAccordingToTheirMarginalProb(
			const int pos);
	doubleRep getNodeHighestMarginal(	const tree::nodeP& inNodeP);
};


#endif
