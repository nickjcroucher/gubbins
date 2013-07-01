#ifndef ___COMPUTE_SUBSTITUTION_COUNTS
#define ___COMPUTE_SUBSTITUTION_COUNTS

#include "definitions.h"
#include "replacementModel.h"
#include "sequenceContainer.h"
#include "tree.h"
#include <map>

class multipleStochasticProcess;
class computeSubstitutionCounts{
public:
	explicit computeSubstitutionCounts(const sequenceContainer& sc, const tree& tr, multipleStochasticProcess* MultSpPtr, string& outDir, VVVdouble& LpostPerSpPerCat, const int simulationsIterNum=1000, const MDOUBLE probCutOffSum=0.5, bool isSilent=false);//DEBUG: Change simulationsIterNum back to 10000

	computeSubstitutionCounts(const computeSubstitutionCounts& other) {*this = other;}	
	computeSubstitutionCounts& operator=(const computeSubstitutionCounts &other);
	virtual ~computeSubstitutionCounts() {}
	void run();
	void computePosteriorOfChangeGivenTerminalsPerSpPerCat();

	void printProbExp();
	void printProbabilityPerPosPerBranch();
	void printProbExpPerPosPerBranch(MDOUBLE probCutOff =0.5,MDOUBLE countsCutOff= 0.2);
	void printExpectationPerBranch();

	void printTreesWithExpectationValuesAsBP(int from,int to);
	void printTreesWithProbabilityValuesAsBP(int from,int to);

	void printProbabilityPerPosPerBranch(int pos, VVVdouble& probChanges, ostream& out, ostream& outCount);
	void printExpectationPerBranch(VVVdouble& expectChanges, ostream& out);
	void printProbExpPerPosPerBranch(int pos, MDOUBLE probCutOff, MDOUBLE countCutOff, VVVdouble& probChanges, VVVdouble& expChanges, ostream& out, ostream& outCount);


	map<int,map<int,vector<double> > > get_expMap_father2son() {return _expMap_father2son;};
	map<int,map<int,vector<double> > > get_probMap_father2son() {return _probMap_father2son;};
	
	VVVVdouble getExpChanges(){return _expChanges_PosNodeXY;};		// expChanges_PosNodeXY[pos][nodeID][x][y]
	VVVVdouble getProbChanges(){return _probChanges_PosNodeXY;};	// probChangesForBranch[pos][nodeID][x][y]
	VVVVdouble getJointProb(){return _jointProb_PosNodeXY;};		// _jointProb_PosNodeXY[pos][nodeID][x][y]	


protected:
//members
	int _alphabetSize;
	const tree _tr;
	const sequenceContainer _sc;

	multipleStochasticProcess* _pMSp;  

	sequence* _refSeq; // the reference sequence
	string _outDir;
	bool _isSilent;
	int _simulationsIterNum;
	MDOUBLE _probCutOffSum;

	VVdouble _LpostPerCat; // the posterior probability for each position for each rate category
	VVVdouble _LpostPerSpPerCat; // _LpostPerSpPerCat[sp][rateCat][pos]


	map<int,map<int,vector<double> > > _expMap_father2son;

	map<int,map<int,vector<double> > > _probMap_father2son;

	//VVVVdouble _posteriorsGivenTerminals;	// posteriorsGivenTerminals[pos][nodeID][x][y]
	VVVVdouble _probChanges_PosNodeXY;		// probChanges_PosNodeXY[pos][nodeID][fatherState][sonState] - after simulations
	VVVVdouble _expChanges_PosNodeXY;		// expChanges_PosNodeXY[pos][nodeID][fatherState][sonState] - after simulations and postProb
	VVVVdouble _jointProb_PosNodeXY;		// probJoint_PosNodeXY[pos][nodeID][fatherState][sonState] - after computePosteriorOfChangeGivenTerminals

};

#endif
