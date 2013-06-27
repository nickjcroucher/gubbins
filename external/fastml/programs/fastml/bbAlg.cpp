#include "bbAlg.h"
#include "computeUpAlg.h"
#include "likelihoodComputation.h"
#include "maseFormat.h"
#include <cmath>

bbAlg::bbAlg(const tree& et,
					vector<stochasticProcess> &spVec,
					const sequenceContainer& sc,
					const bbAlg::boundMethod boundType,
					const string& reportFileName,
					const MDOUBLE computeAgainExactTreshold,
					const distribution * forceDistr) :
	_reportFileName(reportFileName),
	BandBReportAllPos1(reportFileName,et.getInternalNodesNum()*spVec[0].alphabetSize()*sc.seqLen()),
	_et(et), _spVec(spVec), _sc(sc)
{
	cout<<"in bbAlg"<<endl;
	_boundMethod = boundType; 
	_alphabetSize=_spVec[0].alphabetSize();	
	_seqLen=_sc.seqLen();
	if (_spVec.size()>1) {//w codon model + gamma special case
		_cpij._V.resize(forceDistr->categories());
		for (int i=0; i < _spVec.size(); ++i) 
			_cpij._V[i].fillPij(_et,_spVec[i]);	
		_spVec[0].setDistribution(forceDistr);//update the first process with gamma distr 
												//for all the functions that needs number catregor and categor probabilty
	}
	else{
		cout<<"no codon model"<<endl;
		_cpij.fillPij(_et,_spVec[0]); 
	}

	_bbesavp1 = new bbEvaluateSpecificAV(_et,_spVec[0],_sc,_cpij);
					
	_bbNodeOrderAlg1 = new bbNodeOrderAlg(_et,_spVec[0],_sc,_cpij,computeAgainExactTreshold);
	cout<<"after bbNodeOrderAlg"<<endl;
	_bbfindBestAVDynProg1 = new bbfindBestAVDynProg(&_et,&_spVec[0],_sc,&_cpij);
	cout<<"after bbfindBestAVDynProg"<<endl;
	sequence tmp(_sc.getAlphabet());
	const int startingVal = -2;
	tmp.resize(_seqLen,&startingVal);
	cout<<"after resize"<<endl;
	_internalSequences.resize(_et.getNodesNum(),tmp);
	cout<<"after _internalSequences resize"<<endl;
	_bestReconstruction.resize(_et.getNodesNum(),tmp);
	cout<<"afetr _bestReconstruction resize"<<endl;

}

void bbAlg::outputTheJointProbAtEachSite(const string & outputFileProbJoint) {
	ofstream jointProbOutput(outputFileProbJoint.c_str());
	MDOUBLE totalLogLikelihood =0;
	for (int j=0; j < _jointL.size(); ++j) {
		totalLogLikelihood+=log(_jointL[j]);
		jointProbOutput<<"Joint log likelihood of position "<<j+1;// j+1 so that positions start from 1, and not from 0.
		jointProbOutput<<": "<<log(_jointL[j])<<endl;
	}
	jointProbOutput<<"total log likelihood of joint reconstruction: "<<totalLogLikelihood<<endl;
	jointProbOutput.close();
}

MDOUBLE bbAlg::bbReconstructAllPositions(sequenceContainer& res){
	cout<<"in bbAlg::bbReconstructAllPositions"<<endl;
	MDOUBLE sumLogLikelihood=0;
	computePijGam cpij;
	cout<<"Gamma model. Branch and Bound.\nReconstructing position: ";
	_jointL.clear();
	for (int i=0 ; i < _seqLen ; ++i) {
		fillProbOfPosition(i);		
		_bbReport = new BandBReport(_reportFileName,i,_spVec[0].alphabetSize());
		MDOUBLE tmp = bbReconstructPositions(i);
		_jointL.push_back(tmp);
		assert(tmp>0);
		sumLogLikelihood+=log(tmp);
		if (_reportFileName!="") {
			if (_bbReport->size()>20*_et.getInternalNodesNum()) {
				_bbReport->makeReport();
			} else 	if (_bbReport->size()<20*_et.getInternalNodesNum()) {
				errorMsg::reportError("error in function bbReconstructAllPositions");
			}
			BandBReportAllPos1.totalNumberOfNodeVisited += _bbReport->size();
		}
		delete _bbReport;
	}
	res = fromAncestralSequenceToSeqData(); // returning the ancestral sequences
	BandBReportAllPos1.printReport();
	return sumLogLikelihood;
}

MDOUBLE bbAlg::bbReconstructPositions(const int pos){
	_bestRecord=0;
	return bbReconstructPositions(pos,1); // 1 - start the first node in the search tree.

}

MDOUBLE bbAlg::bbReconstructPositions(const int pos,
						 const int nodeNum) {
	tree::nodeP node2check=NULL;
	vector<int> charOrder;
	doubleRep exactVal=0;
	if (nodeNum == 1) {
		_bbNodeOrderAlg1->getNextNodeAndCharOrder(	node2check,
													charOrder,
													_internalSequences,
													pos,
													true,
													exactVal);
	}
	else {
		_bbNodeOrderAlg1->getNextNodeAndCharOrder(	node2check,
													charOrder,
													_internalSequences,
													pos,
													false,
													exactVal);
	}
	int k;
	for (k = 0; k < charOrder.size(); k++) {
		_internalSequences[node2check->id()][pos] = charOrder[k];
		bool haveToGoDown=false;
		if (nodeNum<_et.getInternalNodesNum()) {
			MDOUBLE boundSigma,boundMax;
			haveToGoDown =decideIfHaveToGoDown(pos,boundSigma,boundMax);
        	_bbReport->report(	node2check->name(),
								charOrder[k],
								nodeNum,
								_bestRecord/_pOfPos,
								0.00,
								boundSigma/_pOfPos,
								boundMax/_pOfPos);
		};
		if (haveToGoDown == true) {
			bbReconstructPositions(pos,(nodeNum+1));
		}
	

		if (nodeNum==_et.getInternalNodesNum()) {
			MDOUBLE tmp = _bbesavp1->evaluateSpecificAv(pos,&_internalSequences);
			if (tmp > _bestRecord) {
				vector<tree::nodeP> allNodes;
				_et.getAllHTUs(allNodes,_et.getRoot());
				for (int j = 0 ; j < allNodes.size(); j++) {
					_bestReconstruction[allNodes[j]->id()][pos]=_internalSequences[allNodes[j]->id()][pos];
				}
				_bestRecord = tmp;
			}
			_bbReport->report(	node2check->name(),
								charOrder[k],
								nodeNum,
								_bestRecord/_pOfPos,
								tmp/_pOfPos,
								0.0,
								0.0);
		}
	}

	_internalSequences[node2check->id()][pos] = -2;
	_bbNodeOrderAlg1->putBack(node2check,exactVal);
	return _bestRecord;
}



bbAlg::~bbAlg() {	delete _bbNodeOrderAlg1;
					delete _bbesavp1;
					delete _bbfindBestAVDynProg1;}

void bbAlg::fillProbOfPosition(const int pos) {
	
	_pOfPos = likelihoodComputation::getLofPos(pos,_et,_sc,_cpij,_spVec[0]);
}



sequenceContainer bbAlg::fromAncestralSequenceToSeqData() {
	int j=0;
	sequenceContainer sD;
	for (j=0; j < _sc.numberOfSeqs(); ++j) {
		sD.add(_sc[j]);
	}
	vector<tree::nodeP> HTUs;
	_et.getAllHTUs(HTUs,_et.getRoot());
	for (j=0; j < HTUs.size(); ++j) {
		sequence tmpSeq(_sc.getAlphabet());
		for (int pos=0; pos<_seqLen;++pos) {
			tmpSeq.push_back(_bestReconstruction[HTUs[j]->id()][pos]);
		}
		tmpSeq.setID(sD.numberOfSeqs());
		tmpSeq.setName(HTUs[j]->name());
		sD.add(tmpSeq);
	}
	return sD;
}





bool bbAlg::decideIfHaveToGoDown(const int pos,
								 MDOUBLE& boundSigma,
								 MDOUBLE& boundMax) const {
//---------------------------------------------------------------------
// checkBoundSigma and checkBoundMax return true, if we have to go down
// in the search tree. This is also the ouput of this function.
// i.e., the bound is always an upper bound on the results.
// it is compared with the best score so far, i.e., the lower bound,
// and if the upperbound<lowerbound that there is no need going down.
// When the two bounds are used, 
// it is enough that one is false to indicate no need to go down.
//---------------------------------------------------------------------

	bool acor1 = false;
	bool acor2 = false;
	switch (_boundMethod) {
	case max: return checkBoundMax(pos,boundMax);
					break;
	case sum: return checkBoundSigma(pos,boundSigma);
				break;
	case both: 
			acor1 = checkBoundSigma(pos,boundSigma);
			acor2 = checkBoundMax(pos,boundMax);

//			if ((acor1 == true) && (acor2 == false)) {
//				cerr<<"max is better"<<endl;
//			} else if ((acor2 == true) && (acor1 == false)) {
//				cerr<<"sum is better"<<endl;
//			}  
			return (acor1 && acor2); 
			break;
		default: errorMsg::reportError("Error in function decideIfHaveToGoDown");
	}

	errorMsg::reportError("Error in function decideIfHaveToGoDown");
	return true;
}

bool bbAlg::checkBoundSigma(const int pos,
							MDOUBLE& inBoundSigma) const {
	inBoundSigma = _bbesavp1->evaluateSpecificAv(pos,&_internalSequences);
	if (inBoundSigma < _bestRecord) return false;
	else return true;
}

bool bbAlg::checkBoundMax(const int pos, MDOUBLE& inboundMax) const {
	// to make
	inboundMax = 0.0;
//	MDOUBLE rate;
	for (int rateCategor=0; rateCategor < _spVec[0].categories(); rateCategor++) {
		inboundMax+= (
			_bbfindBestAVDynProg1->evaluateSpecificAvDP(pos,&_internalSequences,rateCategor)*
			_spVec[0].ratesProb(rateCategor));
	}
	if (inboundMax < _bestRecord) return false;
	else return true;
}


