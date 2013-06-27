#include "bbEvaluateSpecificAV.h"

bbEvaluateSpecificAV::bbEvaluateSpecificAV(const tree& et,
					const stochasticProcess& sp,
					const sequenceContainer& sc,
					const computePijGam& cpij) : _et(et), _sp(sp), _sc(sc), _bbcpij(cpij) {
	_sctm = new seqContainerTreeMap(_sc,_et);

	_alphabetSize=_sc.alphabetSize();
	_Lvec.resize(_et.getNodesNum());
	for (int i=0; i < _Lvec.size(); ++i ) {
        _Lvec[i].resize(_alphabetSize);
	}
}

bbEvaluateSpecificAV::~bbEvaluateSpecificAV() {
	delete _sctm;
}

MDOUBLE bbEvaluateSpecificAV::evaluateSpecificAv(
		const int pos,
		const vector<sequence>* ancestralSequences) {
	_ancss = ancestralSequences;
	return recursiveEvaluateSpecificAv(pos,_et.getRoot());
}

MDOUBLE bbEvaluateSpecificAV::recursiveEvaluateSpecificAv(
					const int pos,
					const tree::nodeP thisNode) {

	MDOUBLE res=0.0;
	for (int rateCategor=0;rateCategor<_sp.categories();rateCategor++) {
		res += (
			recursiveEvaluateSpecificAv(pos,thisNode,rateCategor)*
			_sp.ratesProb(rateCategor)
			);
	}
	return res;
}

MDOUBLE bbEvaluateSpecificAV::recursiveEvaluateSpecificAv(const int pos,
											const tree::nodeP thisNode,
											const int categor) {
	
	int letterInNode;
	if (thisNode->isLeaf() ) {
		const int seqID = _sctm->seqIdOfNodeI(thisNode->id());
		letterInNode = _sc[seqID][pos];
		for (int k = 0; k < _alphabetSize ; ++k) { // taking care of ? by the -2 64 - for codons...
			if ((letterInNode==-2) || (letterInNode==-1)||(letterInNode==64) ||(letterInNode==k)) _Lvec[thisNode->id()][k] = 1.0;
			else _Lvec[thisNode->id()][k] = 0.0;
		}
		return 0.0;
	}

	for (int i = 0 ; i < thisNode->getNumberOfSons() ; ++i ) {// recursive call for the childs
		recursiveEvaluateSpecificAv(pos,thisNode->getSon(i),categor);
	}

	letterInNode = (*_ancss)[thisNode->id()][pos];
	if (letterInNode == -2) {// internal node with asterix.
		for (int y = 0 ; y < _alphabetSize ; ++y) {
			MDOUBLE rate = _sp.rates(categor); // the r.
			_Lvec[thisNode->id()][y] = 1.0;
			for (int u = 0 ; u < thisNode->getNumberOfSons() ; ++u) {
				MDOUBLE tmp = 0;
				for (int letInSon = 0 ; letInSon<_alphabetSize; ++letInSon) {
					tmp+=(
						_bbcpij.getPij(categor,thisNode->getSon(u)->id(),y,letInSon)*
						_Lvec[thisNode->getSon(u)->id()][letInSon]
						);
				}
				_Lvec[thisNode->id()][y] *= tmp;
				
			}
		}
	}
	
	else { // if the character in the HTU is known (not an asterix)
		for (int w = 0 ; w < _alphabetSize ; ++w) {
			if (w != letterInNode) _Lvec[thisNode->id()][w] = 0.0;
			else {
//				MDOUBLE rate = _myStoc_proc.rates(categor); // the r.
				_Lvec[thisNode->id()][w] = 1.0;
				for (int z = 0 ; z < thisNode->getNumberOfSons() ; ++z) {
					MDOUBLE tmp = 0;
					for (int letInSon = 0 ; letInSon<_alphabetSize; ++letInSon) {
						tmp += (
							_bbcpij.getPij(categor,thisNode->getSon(z)->id(),w,letInSon)*
							_Lvec[thisNode->getSon(z)->id()][letInSon]
							);
					}
					_Lvec[thisNode->id()][w] *= tmp;
				}
			}// end of else
		}
	}
	
	MDOUBLE result= 0.0;
	if (thisNode->father() ==  NULL){ // tree root
		
		for (int letRoot = 0 ; letRoot < _alphabetSize; ++letRoot) {
			result += _sp.freq(letRoot) * _Lvec[thisNode->id()][letRoot];
		}
	}
	return result;

}





