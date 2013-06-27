#include "sequenceDataDiff.h"
#include <iostream>
using namespace std;

void sequenceDataDiff::computeDifferences(){
	for (int i=0;i<_sc1.numberOfSeqs();++i) {
		string name1 = _sc1[i].name();
		int idOf1in2 = _sc2.getId(name1,false);//return -1 if not found...
		if (idOf1in2==-1) {
			string x = "sequence does not exist ";
			x+=name1;
			unitDiff ud(x);
			_differences.push_back(ud);
			continue;
		}
		const sequence& sequence1 = _sc1[i];
		const sequence& sequence2 = _sc2[i];
		if (sequence1.seqLen() != sequence1.seqLen()) {
			string x = "sequences don't have the same length ";
			x+=name1;
			unitDiff ud(x);
			_differences.push_back(ud);
			continue;
		}

		for (int j=0; j < sequence1.seqLen(); ++j) {
			if (sequence1[j] != sequence2[j]) {
				unitDiff ud(name1,j,sequence1.toString(j),sequence2.toString(j));
				_differences.push_back(ud);
			}
		}
	}
}


void sequenceDataDiff::printDiff(ostream& out) {
	for (int i=0; i < _differences.size(); ++i) {
		out<<_differences[i]._seqName;
		out<<" ";
		out<<_differences[i]._pos;
		out<<" ";
		out<<_differences[i]._letInSd1;
		out<<" ";
		out<<_differences[i]._letInSd2;
		out<<endl;
	}
}


