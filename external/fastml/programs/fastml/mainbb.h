#ifndef ___BB__MAIN__FILE
#define ___BB__MAIN__FILE

#include "bb_options.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "tree.h"
#include "codon.h"
#include "nucleotide.h"

#include "suffStatComponent.h"

#include <vector>
using namespace std;


class mainbb {
public:
	explicit mainbb(int argc, char* argv[]);
	virtual ~mainbb();

private:
	const bb_options* _options;
	sequenceContainer _sc;
	sequenceContainer _originSc; //hold the sc before change the gaps
	tree _et;
	vector<stochasticProcess> _spVec; //hold stochastic process 
									//if codon yang model with gamma then 
									//holds number of categores of replacment model 
	distribution *_forceDistr; //holds the w distribution of yang codon model. 

	alphabet* _alph;
	sequenceContainer _resulutingJointReconstruction;

	void getStartingStochasticProcess();
	void createStochasticProcessVec();
	Vdouble computeFreq(codon &codonAlph);
	Vdouble computeGTRFreq(nucleotide &nucAlph);
	Vdouble freqGTR(const sequenceContainer &nucSc, nucleotide * nucAlph);
	
	// get starting tree
	void getStartingEvolTreeTopology();
	void getStartingNJtreeNjMLdis();
	void getStartingTreeNJ_fromDistances(const VVdouble& disTab,const vector<string>& vNames);
	void getStartingTreeFromTreeFile();
	void getStartingBranchLengthsAndAlpha();
	void printOutputTree();

	//get starting tree and codon model
	 void getStartingBLAndModelParam();

	// JOINT WITH GAMMA
	void printAncestralSequencesGammaJoint();
	void findAncestralSequencesGammaJoint();

	// JOINT WITHOUT GAMMA
	void findAncestralSequencesHomJoint();

	// MARGINAL RECONSTRUCTION:
	void getMarginalReconstruction();


	void fillOptionsParameters(int argc, char* argv[]);
	void getStartingSequenceData();
	void printSearchParameters();
	void printBBProjectInfo();
	void replaceSequences(sequenceContainer &sc2change,sequenceContainer &originSc);


};


#endif

