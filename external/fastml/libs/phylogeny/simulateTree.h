// $Id: simulateTree.h 3574 2008-02-27 10:44:30Z itaymay $

#ifndef ___SIMULATE_TREE
#define ___SIMULATE_TREE

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"

//class sequenceData; // to be able to go to simulate data.

class simulateTree {
public:
	explicit simulateTree(const tree&  _inEt,const stochasticProcess& sp,
		const alphabet* alph);
	void generate_seq(int seqLength);
	
	// This function generates the sequences not using the discrete gamma, but rather,
	// the rates are sampled from the continuous distribution.
	// It assumes the Gamma distribution has mean 1 (alpha = beta).
	void generate_seq_continuous_gamma(int seqLength);

	void generate_seqWithRateVector(const Vdouble& simRates, const int seqLength);	
	//these function do the same simulation as above but check that no stop codon is created.
	//applicable only when the stochasticProcess is based on nucleotides
	void generate_seqWithRateVectorNoStopCodon(const Vdouble& simRates, int seqLength);

	tree gettree() {return _et;}
	virtual ~simulateTree();
	sequenceContainer toSeqData();
	sequenceContainer toSeqDataWithoutInternalNodes();
	void generate_rates_continuous_gamma(const int seqLength,const MDOUBLE alpha,Vdouble rates);

private:
	void generateRootSeq(int seqLength);	
	void recursiveGenerateSpecificSeq(const Vdouble& rateVec, int seqLength, tree::nodeP myNode);
	int giveRandomChar() const;
	int giveRandomChar(const int letterInFatherNode, const MDOUBLE length,const int pos) const;
	int getRandCategory(const int pos) const;

	vector<sequence> _simulatedSequences; // the sequences (nodes * seqLen)
	tree _et;
	const stochasticProcess& _sp;
	const alphabet* _alph;
};

#endif

