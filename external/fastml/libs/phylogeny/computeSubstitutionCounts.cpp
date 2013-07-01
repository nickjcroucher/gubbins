#include "computeSubstitutionCounts.h"
#include "computePosteriorExpectationOfSubstitutions.h"
#include "computePosteriorExpectationOfSubstitutions_nonReversibleSp.h"
#include "multipleStochasticProcess.h"
#include "matrixUtils.h"
#include "simulateJumps.h"
#include "simulateCodonsJumps.h"
#include "simulateJumpsAbstract.h"
#include "treeIt.h"
#include "treeUtil.h"

/********************************************************************************************
computeSubstitutionCounts
*********************************************************************************************/
computeSubstitutionCounts::computeSubstitutionCounts(const sequenceContainer& sc, const tree& tr, multipleStochasticProcess* MultSpPtr, string& outDir, VVVdouble& LpostPerSpPerCat, const int simulationsIterNum, const MDOUBLE probCutOffSum, bool isSilent):
_tr(tr),_sc(sc),_pMSp(MultSpPtr),_outDir(outDir),_LpostPerSpPerCat(LpostPerSpPerCat), _simulationsIterNum(simulationsIterNum), _probCutOffSum(probCutOffSum),_isSilent(isSilent)
{
	if(!_pMSp->getSPVecSize()){
		errorMsg::reportError("Trying to call computeSubstitutionCounts with an empty multipleStochasticProcess object at computeSubstitutionCounts::computeSubstitutionCounts");
	}
	_alphabetSize = _pMSp->getSp(0)->alphabetSize();
}

computeSubstitutionCounts& computeSubstitutionCounts::operator=(const computeSubstitutionCounts &other){
	if (this != &other) {              // Check for self-assignment
	}
	return *this;
}


/********************************************************************************************
*********************************************************************************************/
void computeSubstitutionCounts::run()
{
	for(int fatherStateIndex = 0;fatherStateIndex < _alphabetSize;++fatherStateIndex){
		for(int sonStateIndex = 0;sonStateIndex < _alphabetSize;++sonStateIndex){
			//if(sonStateIndex == fatherStateIndex) continue;
			_expMap_father2son[fatherStateIndex][sonStateIndex].resize(_sc.seqLen(),0);
			_probMap_father2son[fatherStateIndex][sonStateIndex].resize(_sc.seqLen(),0);
		}
	}

	resize_VVVV(_sc.seqLen(),_tr.getNodesNum(),_alphabetSize,_alphabetSize,_jointProb_PosNodeXY);
	resize_VVVV(_sc.seqLen(),_tr.getNodesNum(),_alphabetSize,_alphabetSize,_probChanges_PosNodeXY);
	resize_VVVV(_sc.seqLen(),_tr.getNodesNum(),_alphabetSize,_alphabetSize,_expChanges_PosNodeXY);

	computePosteriorOfChangeGivenTerminalsPerSpPerCat();	// GLM - multiple SPs
}

/********************************************************************************************
*********************************************************************************************/
void computeSubstitutionCounts::computePosteriorOfChangeGivenTerminalsPerSpPerCat()
{	
	int numOfSPs = _pMSp->getSPVecSize();

	// per Sp
	for (int spIndex=0; spIndex < numOfSPs; ++spIndex) {
		// Per RateCategory -- All the computations are done while looping over rate categories
		stochasticProcess * currentSp = _pMSp->getSp(spIndex);
		int numOfRateCategories = currentSp->categories();	
		for (int rateCategIndex=0 ; rateCategIndex < numOfRateCategories;++rateCategIndex)
		{
			tree copy_et = _tr;
			MDOUBLE rateCategVal = currentSp->rates(rateCategIndex);
			MDOUBLE  minimumRateCategVal = 0.0000001;
			MDOUBLE rate2multiply = max(rateCategVal,minimumRateCategVal);
			if(rateCategVal < minimumRateCategVal){
				LOGnOUT(4, <<" >>> NOTE: the rate category "<<rateCategVal<<" is too low for computePosteriorExpectationOfChangePerSite"<<endl);	}
			copy_et.multipleAllBranchesByFactor(rate2multiply);
			//if(!_isSilent) 
				//LOGnOUT(4, <<"running "<<gainLossOptions::_numOfSimulationsForPotExp<<" simulations for rate "<<rate2multiply<<endl);
			simulateJumpsAbstract* simPerRateCategory;
			if(_alphabetSize == 61)
				simPerRateCategory = new simulateCodonsJumps(copy_et,*currentSp,_alphabetSize);
			else
				simPerRateCategory = new simulateJumps(copy_et,*currentSp,_alphabetSize);
				
			simPerRateCategory->runSimulation(_simulationsIterNum);
			if(!_isSilent) 
				LOGnOUT(4,<<"finished simulations"<<endl);

			// Per POS		
			for (int pos = 0; pos <_sc.seqLen(); ++pos)
			{
				LOG(6,<<"pos "<<pos+1<<endl);
				// I) computePosteriorOfChangeGivenTerminals
				VVVdouble posteriorsGivenTerminalsPerRateCategoryPerPos;
				computePosteriorExpectationOfSubstitutions* cpesPerRateCategoryPerPos ;
				if(currentSp->isReversible())
					cpesPerRateCategoryPerPos = new computePosteriorExpectationOfSubstitutions(copy_et,_sc,currentSp);	// Per POS,CAT
				else
					cpesPerRateCategoryPerPos = new computePosteriorExpectationOfSubstitutions_nonReversibleSp(copy_et,_sc,currentSp);	// Per POS,CAT
				cpesPerRateCategoryPerPos->computePosteriorOfChangeGivenTerminals(posteriorsGivenTerminalsPerRateCategoryPerPos,pos);

				// II) Exp - take in account both: 1) simulations 2) posteriorsGivenTerminal
				VVVdouble expChangesForBranchPerRateCategoryPerPos;	// Sim+Exp
				resize_VVV(_tr.getNodesNum(),_alphabetSize,_alphabetSize,expChangesForBranchPerRateCategoryPerPos);

				VVdouble expVV = cpesPerRateCategoryPerPos->computeExpectationAcrossTree(*simPerRateCategory,posteriorsGivenTerminalsPerRateCategoryPerPos,
					expChangesForBranchPerRateCategoryPerPos);	// Per POS
				for(int fatherStateIndex = 0;fatherStateIndex < _alphabetSize;++fatherStateIndex){
					for(int sonStateIndex = 0;sonStateIndex < _alphabetSize;++sonStateIndex){
						if(sonStateIndex == fatherStateIndex) continue;
						_expMap_father2son[fatherStateIndex][sonStateIndex][pos] += expVV[fatherStateIndex][sonStateIndex]*_LpostPerSpPerCat[spIndex][rateCategIndex][pos];
					}
				}

				// III) Sim - take in account both: 1) simulations 2) posteriorsGivenTerminal
				VVVdouble probChangesForBranchPerRateCategoryPerPos;	// Sim+Prob
				resize_VVV(_tr.getNodesNum(),_alphabetSize,_alphabetSize,probChangesForBranchPerRateCategoryPerPos);
				VVdouble probVV = cpesPerRateCategoryPerPos->computePosteriorAcrossTree(*simPerRateCategory,posteriorsGivenTerminalsPerRateCategoryPerPos,probChangesForBranchPerRateCategoryPerPos);
				for(int fatherStateIndex = 0;fatherStateIndex < _alphabetSize;++fatherStateIndex){
					for(int sonStateIndex = 0;sonStateIndex < _alphabetSize;++sonStateIndex){
						if(sonStateIndex == fatherStateIndex) continue;
						_probMap_father2son[fatherStateIndex][sonStateIndex][pos] += probVV[fatherStateIndex][sonStateIndex]*_LpostPerSpPerCat[spIndex][rateCategIndex][pos];
					}
				}
				//	Store all information PerCat,PerPOS
				for(int i=0;i<_probChanges_PosNodeXY[pos].size();++i){	// nodeId
					for(int j=0;j<_probChanges_PosNodeXY[pos][i].size();++j){	// fatherState
						for(int k=0;k<_probChanges_PosNodeXY[pos][i][j].size();++k){	// sonState
							_jointProb_PosNodeXY[pos][i][j][k] += posteriorsGivenTerminalsPerRateCategoryPerPos[i][j][k]*_LpostPerSpPerCat[spIndex][rateCategIndex][pos];
							_probChanges_PosNodeXY[pos][i][j][k] += probChangesForBranchPerRateCategoryPerPos[i][j][k]*_LpostPerSpPerCat[spIndex][rateCategIndex][pos];
							_expChanges_PosNodeXY[pos][i][j][k] += expChangesForBranchPerRateCategoryPerPos[i][j][k]*_LpostPerSpPerCat[spIndex][rateCategIndex][pos];
						}
					}
				}
				delete(cpesPerRateCategoryPerPos);
			}
			delete(simPerRateCategory);
			// Per POS
		}
		// per rateCat
	}
	// Per Sp
}



/********************************************************************************************
printProbExp()
print perPos (over all branches)
use the members _expV01, _expV10 for basic 
*********************************************************************************************/
void computeSubstitutionCounts::printProbExp()
{

	string posteriorExpectationOfChangeString = _outDir + "//" + "posteriorExpectationOfChange.txt";
	ofstream posteriorExpectationStream(posteriorExpectationOfChangeString.c_str());
	string posteriorProbabilityOfChangeString = _outDir + "//" + "posteriorProbabilityOfChange.txt";
	ofstream posteriorProbabilityStream(posteriorProbabilityOfChangeString.c_str());

	int fatherStateIndex,sonStateIndex;
	posteriorExpectationStream<<"#POS"<<"\t";
	posteriorProbabilityStream<<"#POS"<<"\t";

	for (fatherStateIndex = 0;fatherStateIndex < _alphabetSize;++fatherStateIndex){
		for (sonStateIndex = 0;sonStateIndex < _alphabetSize;++sonStateIndex){
			if(sonStateIndex == fatherStateIndex) continue;
			posteriorExpectationStream<<_sc.getAlphabet()->fromInt(fatherStateIndex)<<"->"<<_sc.getAlphabet()->fromInt(sonStateIndex)<<"\t";
			posteriorProbabilityStream<<_sc.getAlphabet()->fromInt(fatherStateIndex)<<"->"<<_sc.getAlphabet()->fromInt(sonStateIndex)<<"\t";
		}
	}
	posteriorExpectationStream<<endl;
	posteriorProbabilityStream<<endl;

	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		posteriorExpectationStream<<pos+1<<"\t";
		posteriorProbabilityStream<<pos+1<<"\t";
		for (fatherStateIndex = 0;fatherStateIndex < _alphabetSize;++fatherStateIndex){
			for (sonStateIndex = 0;sonStateIndex < _alphabetSize;++sonStateIndex){
				if(sonStateIndex == fatherStateIndex) continue;//ofir, note the change in print format
				posteriorExpectationStream<<_expMap_father2son[fatherStateIndex][sonStateIndex][pos]<<"\t";
				posteriorProbabilityStream<<_probMap_father2son[fatherStateIndex][sonStateIndex][pos]<<"\t";	
			}
		}
		posteriorExpectationStream<<endl;
		posteriorProbabilityStream<<endl;
	}
	posteriorExpectationStream.close();
	posteriorProbabilityStream.close();
}


/********************************************************************************************
printProbabilityPerPosPerBranch 1
produce 2 print files:
1. print detailed file (out)
2. print summary over all branches (outSum)
*********************************************************************************************/
void computeSubstitutionCounts::printProbabilityPerPosPerBranch()
{
	string probabilityPerPosPerBranch = _outDir + "//" + "probabilityPerPosPerBranch.txt"; 
	ofstream probabilityPerPosPerBranchStream(probabilityPerPosPerBranch.c_str());
	probabilityPerPosPerBranchStream<<"# print values over probCutOff "<<_probCutOffSum<<endl;
	probabilityPerPosPerBranchStream<<"#Event"<<"\t"<<"POS"<<"\t"<<"branch"<<"\t"<<"branchLength"<<"\t"<<"distance2root"<<"\t"<<"probability"<<endl;
	
	string countProbPerPos = _outDir + "//" + "probabilityPerPos.txt"; 
	ofstream countProbPerPosStream(countProbPerPos.c_str());
	countProbPerPosStream<<"# print values over probCutOff "<<_probCutOffSum<<endl;
	countProbPerPosStream<<"#POS"<<"\t";
	for(int fatherStateIndex = 0;fatherStateIndex < _alphabetSize;++fatherStateIndex){
		for(int sonStateIndex = 0;sonStateIndex < _alphabetSize;++sonStateIndex){
			if(sonStateIndex == fatherStateIndex) continue;
			countProbPerPosStream<<"prob"<<_sc.getAlphabet()->fromInt(fatherStateIndex)<<"->"<<_sc.getAlphabet()->fromInt(sonStateIndex)<<"\t";
		}
	}
	countProbPerPosStream<<endl;
	
	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		printProbabilityPerPosPerBranch(pos, _probChanges_PosNodeXY[pos],probabilityPerPosPerBranchStream,countProbPerPosStream);
	}
}
/********************************************************************************************
printGainLossProbabilityPerPosPerBranch 1.1
*********************************************************************************************/
void computeSubstitutionCounts::printProbabilityPerPosPerBranch(int pos, VVVdouble& probChanges, ostream& out, ostream& outCount)
{
	VVdouble countFromFather2Son;
	countFromFather2Son.resize(_alphabetSize);
	int fatherStateIndex,sonStateIndex;
	treeIterTopDownConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		for(fatherStateIndex = 0;fatherStateIndex < _alphabetSize;++fatherStateIndex){
			countFromFather2Son[fatherStateIndex].resize(_alphabetSize,0);
			for(sonStateIndex = 0;sonStateIndex < _alphabetSize;++sonStateIndex){
				if(sonStateIndex == fatherStateIndex) continue;
				if(probChanges[mynode->id()][fatherStateIndex][sonStateIndex] > _probCutOffSum){//NIM
					out<<_sc.getAlphabet()->fromInt(fatherStateIndex)<<"->"<<_sc.getAlphabet()->fromInt(sonStateIndex)<<"\t"<<pos+1<<"\t"<<mynode->name()<<"\t"<<mynode->dis2father()<<"\t"<<getDistanceFromNode2ROOT(mynode)<<"\t"<<probChanges[mynode->id()][fatherStateIndex][sonStateIndex]<<endl;
					countFromFather2Son[fatherStateIndex][sonStateIndex] += probChanges[mynode->id()][fatherStateIndex][sonStateIndex];
				}
			}
		}
	}
	outCount<<pos+1<<"\t";
	for(fatherStateIndex = 0;fatherStateIndex < _alphabetSize;++fatherStateIndex){
		for(sonStateIndex = 0;sonStateIndex < _alphabetSize;++sonStateIndex){
			if(sonStateIndex == fatherStateIndex) continue;
			//if(countFromFather2Son[fatherStateIndex][sonStateIndex] == 0) continue;//NIMROD
			outCount<<countFromFather2Son[fatherStateIndex][sonStateIndex]<<"\t";
		}
	}
	outCount<<endl;
}



/********************************************************************************************
*********************************************************************************************/
void computeSubstitutionCounts::printExpectationPerBranch()
{
	// ExpectationPerBranch
	VVVdouble posteriorsGivenTerminalsTotal;
	resize_VVV(_tr.getNodesNum(),_alphabetSize,_alphabetSize,posteriorsGivenTerminalsTotal);
	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		for(int i=0;i<_expChanges_PosNodeXY[pos].size();++i){
			for(int j=0;j<_expChanges_PosNodeXY[pos][i].size();++j){
				for(int k=0;k<_expChanges_PosNodeXY[pos][i][j].size();++k){
					posteriorsGivenTerminalsTotal[i][j][k] += _expChanges_PosNodeXY[pos][i][j][k];
				}
			}
		}
	}
	string expectationPerBranch = _outDir + "//" + "ExpectationPerBranch.txt"; 
	ofstream expectationPerBranchStream(expectationPerBranch.c_str());
	printExpectationPerBranch(posteriorsGivenTerminalsTotal,expectationPerBranchStream);
}
/********************************************************************************************
*********************************************************************************************/
void computeSubstitutionCounts::printExpectationPerBranch(VVVdouble& expectChanges, ostream& out)
{
	treeIterTopDownConst tIt(_tr);
	out<<"#Event"<<"\t"<<"branch"<<"\t"<<"branchLength"<<"\t"<<"distance2root"<<"\t"<<"expectation"<<endl;
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		for(int fatherStateIndex = 0;fatherStateIndex < _alphabetSize;++fatherStateIndex){
			for(int sonStateIndex = 0;sonStateIndex < _alphabetSize;++sonStateIndex){
				if(sonStateIndex == fatherStateIndex) continue;
				out<<_sc.getAlphabet()->fromInt(fatherStateIndex)<<"->"<<_sc.getAlphabet()->fromInt(sonStateIndex)<<"\t"<<
					mynode->name()<<"\t"<<mynode->dis2father()<<"\t"<<getDistanceFromNode2ROOT(mynode)<<"\t"<<expectChanges[mynode->id()][fatherStateIndex][sonStateIndex]<<endl;
			}
		}
	}
}


/********************************************************************************************
*********************************************************************************************/
void computeSubstitutionCounts::printTreesWithExpectationValuesAsBP(int from,int to)
{
	// ExpectationPerPosPerBranch - Print Trees
	Vstring Vnames;
	fillAllNodesNames(Vnames,_tr);
	createDir(_outDir, "TreesWithExpectationValuesAsBP");
	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		string strTreeNum = _outDir + "//" + "TreesWithExpectationValuesAsBP" + "//" + "expTree" + int2string(pos+1) + ".ph";
		ofstream tree_out(strTreeNum.c_str());
		printTreeWithValuesAsBP(tree_out,_tr,Vnames,&_expChanges_PosNodeXY[pos],from,to);
	}
}

/********************************************************************************************
*********************************************************************************************/
void computeSubstitutionCounts::printTreesWithProbabilityValuesAsBP(int from,int to)
{
	// ProbabilityPerPosPerBranch - Print Trees
	Vstring Vnames;
	fillAllNodesNames(Vnames,_tr);
	createDir(_outDir, "TreesWithProbabilityValuesAsBP");
	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		string strTreeNum = _outDir + "//" + "TreesWithProbabilityValuesAsBP"+ "//" + "probTree" + int2string(pos+1) + ".ph";
		ofstream tree_out(strTreeNum.c_str());
		printTreeWithValuesAsBP(tree_out,_tr,Vnames,&_probChanges_PosNodeXY[pos],from,to);
	}
}

/********************************************************************************************
printProbExpPerPosPerBranch 1
produce 2 print files:
1. print detailed file (out)
2. print summary over all branches (outSum)
*********************************************************************************************/
void computeSubstitutionCounts::printProbExpPerPosPerBranch(MDOUBLE probCutOff, MDOUBLE countsCutOff)
{
	string probExpPerPosPerBranch = _outDir + "//" + "expPerPosPerBranch.txt"; 
	ofstream probExpPerPosPerBranchStream(probExpPerPosPerBranch.c_str());
	probExpPerPosPerBranchStream<<"# print values over probCutOff "<<probCutOff<<endl;
	probExpPerPosPerBranchStream<<"#Event"<<"\t"<<"POS"<<"\t"<<"branch"<<"\t"<<"branchLength"<<"\t"<<"distance2root"<<"\t"<<"probability"<<"\t"<<"expectation"<<endl;

	string probExpPerPos = _outDir + "//" + "probExpCountPerPos.txt"; 
	ofstream countProbPerPosStream(probExpPerPos.c_str());
	countProbPerPosStream<<"# print count over probCutOff "<<countsCutOff<<endl;
	countProbPerPosStream<<"#POS"<<"\t"<<"Event"<<"\t"<<"EventProb"<<"\t"<<"EventExp"<<"\t"<<"EventCount"<<endl;
	
	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		printProbExpPerPosPerBranch(pos, probCutOff,countsCutOff, _probChanges_PosNodeXY[pos],_expChanges_PosNodeXY[pos],probExpPerPosPerBranchStream,countProbPerPosStream);
	}
}
/********************************************************************************************
printGainLossProbExpPerPosPerBranch 1.1
Get pos, and iterate over all branches:
1. print detailed file (out)
2. print summary over all branches (outSum)
*********************************************************************************************/
void computeSubstitutionCounts::printProbExpPerPosPerBranch(int pos, MDOUBLE probCutOff, MDOUBLE countCutOff, VVVdouble& probChanges, VVVdouble& expChanges, ostream& out, ostream& outSum)
{
	VVdouble probFather2Son,expFather2Son;
	VVint countFather2Son;
	probFather2Son.resize(_alphabetSize);
	expFather2Son.resize(_alphabetSize);
	countFather2Son.resize(_alphabetSize);
	int fatherStateIndex,sonStateIndex; 

	treeIterTopDownConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		for(fatherStateIndex = 0;fatherStateIndex < _alphabetSize;++fatherStateIndex){
			probFather2Son[fatherStateIndex].resize(_alphabetSize,0);
			expFather2Son[fatherStateIndex].resize(_alphabetSize,0);
			countFather2Son[fatherStateIndex].resize(_alphabetSize,0);
			for(sonStateIndex = 0;sonStateIndex < _alphabetSize;++sonStateIndex){
				if(sonStateIndex == fatherStateIndex) continue;
				out<<_sc.getAlphabet()->fromInt(fatherStateIndex)<<"->"<<_sc.getAlphabet()->fromInt(sonStateIndex)<<"\t"<<
					pos+1<<"\t"<<mynode->name()<<"\t"<<mynode->dis2father()<<"\t"<<getDistanceFromNode2ROOT(mynode)<<"\t"<<probChanges[mynode->id()][fatherStateIndex][sonStateIndex]<<"\t"<<expChanges[mynode->id()][fatherStateIndex][sonStateIndex]<<endl;
				probFather2Son[fatherStateIndex][sonStateIndex] += probChanges[mynode->id()][fatherStateIndex][sonStateIndex];
				expFather2Son[fatherStateIndex][sonStateIndex] += expChanges[mynode->id()][fatherStateIndex][sonStateIndex];
				if (probChanges[mynode->id()][fatherStateIndex][sonStateIndex] > countCutOff)
					countFather2Son[fatherStateIndex][sonStateIndex] += 1;
			}
		}
	}
	for(fatherStateIndex = 0;fatherStateIndex < _alphabetSize;++fatherStateIndex){
		for(sonStateIndex = 0;sonStateIndex < _alphabetSize;++sonStateIndex){
			if(sonStateIndex == fatherStateIndex) continue;
			outSum<<pos+1<<"\t"<<_sc.getAlphabet()->fromInt(fatherStateIndex)<<"->"<<_sc.getAlphabet()->fromInt(sonStateIndex)<<"\t"<<
				probFather2Son[fatherStateIndex][sonStateIndex]<<"\t"<<expFather2Son[fatherStateIndex][sonStateIndex]<<"\t"<<countFather2Son[fatherStateIndex][sonStateIndex]<<endl;
		}
	}
}

