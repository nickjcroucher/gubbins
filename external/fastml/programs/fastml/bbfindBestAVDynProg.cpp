#include "bbfindBestAVDynProg.h"

bbfindBestAVDynProg::bbfindBestAVDynProg(const tree* et,
					const stochasticProcess *sp,
					const sequenceContainer& sc,
					const computePijGam* cpij): _sc(sc) {
	_et = et;
	_sp = sp;
	_bbcpij = cpij;
	_sctm = new seqContainerTreeMap(_sc,*_et);
	_alphabetSize=_sp->alphabetSize();
	_jointLval.resize(_et->getNodesNum());
	_jointCval.resize(_et->getNodesNum());
	for (int i=0; i < _et->getNodesNum(); ++i) {
		_jointLval[i].resize(_alphabetSize);
		_jointCval[i].resize(_alphabetSize);
	}
}

bbfindBestAVDynProg::~bbfindBestAVDynProg() {
	delete _sctm;
}

MDOUBLE bbfindBestAVDynProg::evaluateSpecificAvDP(
		const int pos,
		const vector<sequence>* ancestralSequences,
		const int rateCategor) {
	_ancss = ancestralSequences;

	recursiveComputeLandC(pos,_et->getRoot(),rateCategor);
// modified from NancestralTree::findBestLetInRoot(const int pos) {
	MDOUBLE bestLinRoot =0 ;
	//MDOUBLE bestLetInRoot = -2;
	MDOUBLE tmp = 0.0;
	int letInRoot = (*_ancss)[_et->getRoot()->id()][pos];
	if (letInRoot==-2) {
	
		for (int x = 0 ; x < _alphabetSize; ++x) {
			tmp = _sp->freq(x);
			for (int y =0 ; y < _et->getRoot()->getNumberOfSons() ; ++y) {
				tmp *= _jointLval[_et->getRoot()->getSon(y)->id()][x];
			}
			if (tmp > bestLinRoot) {
				bestLinRoot = tmp;
				//bestLetInRoot = x;
			}
		}
	}
	else {//if (letInRoot!=-2)
		tmp = _sp->freq(letInRoot);
		for (int y =0 ; y < _et->getRoot()->getNumberOfSons() ; ++y) {
			tmp *= _jointLval[_et->getRoot()->getSon(y)->id()][letInRoot];
		}
		if (tmp > bestLinRoot) {
			bestLinRoot = tmp;
			//bestLetInRoot = x;
		}
	}

	//iRoot()->data()[pos] = bestLetInRoot;
	return bestLinRoot;
}

void bbfindBestAVDynProg::recursiveComputeLandC(const int pos,
												const tree::nodeP inNode,
												const int rateCategor) {
// root has to be internal node here.
	for (int i=0; i<inNode->getNumberOfSons();++i) {
		recursiveComputeLandC(pos,inNode->getSon(i),rateCategor);
	}
	if (inNode->father() ==  NULL) return;

	int	letInNode;
	if (inNode->isLeaf()) {
		const int seqID = _sctm->seqIdOfNodeI(inNode->id());
		letInNode=_sc[seqID][pos];
	}
	else {
		letInNode = (*_ancss)[inNode->id()][pos];
	}

	if (letInNode!=-2){ // known leaf, or known HTU, (no root)
		
		for (int FatherLet = 0; FatherLet<_alphabetSize;++FatherLet) {
			_jointLval[inNode->id()][FatherLet] = _bbcpij->getPij(rateCategor,inNode->id(),FatherLet,letInNode);
			_jointCval[inNode->id()][FatherLet] = letInNode;
			for (int k=0; k < inNode->getNumberOfSons() ; ++k) {
				_jointLval[inNode->id()][FatherLet] *= _jointLval[inNode->getSon(k)->id()][letInNode];
			}
		}
	}
	else {// unknown leaf or HTU -> no root.
		for (int letInFather = 0; letInFather < _alphabetSize; ++letInFather) {
			MDOUBLE bestVal = 0;
			int bestLet = -2;
			for (int lenInNode = 0; lenInNode < _alphabetSize; ++lenInNode) {
				MDOUBLE tmp = 1;
				if (inNode->isInternal()) 
					tmp*= _bbcpij->getPij(rateCategor,inNode->id(),letInFather,lenInNode);
				// if it is a leaf, and since it is ? tmp will be 1.0...
				for (int k=0; k < inNode->getNumberOfSons() ; ++k) {
					tmp *= _jointLval[inNode->getSon(k)->id()][lenInNode];
				}
				if (tmp > bestVal) {
					bestVal = tmp;
					bestLet = lenInNode;
				}
			}
			_jointLval[inNode->id()][letInFather] = bestVal;
			_jointCval[inNode->id()][letInFather] = bestLet;
		}
	}
}



