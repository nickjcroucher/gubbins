#ifndef ___SEQ__DATA__DIF
#define ___SEQ__DATA__DIF

#include "sequenceContainer.h"

#include <fstream>
#include <iostream>
#include <string>
using namespace std;

// this class represents a single difference between a pair of sequences.
// I.e., it is used here, to show a difference between two approaches for ancestral sequence
// reconstruction, for example, Joint vs. Marginal, or With and Without Gamma.

class unitDiff{
	friend class sequenceDataDiff;
public:
	explicit unitDiff(const string& seqName,const int pos, const string letInSd1,const string letInSd2) {
		_seqName = seqName; _pos = pos; _letInSd1 = letInSd1; _letInSd2 = letInSd2;
	}
	explicit unitDiff(const string& seqName) { //  in case one seq is only in one
		_seqName = seqName; _pos = -1; _letInSd1 = '?'; _letInSd2 = '?';
	}
private:
	string _seqName;
	int _pos;
	string _letInSd1;
	string _letInSd2;
};

// This class prints differences between two reconstructions (or in general, between any two sequence conatiners)

class sequenceDataDiff {
public:
	sequenceDataDiff(const sequenceContainer& sc1, const sequenceContainer& sc2) :_sc1(sc1) ,_sc2(sc2) {}
	void computeDifferences();
	void printDiff(ostream& out);
private:	
	vector<unitDiff> _differences;
	const sequenceContainer& _sc1;
	const sequenceContainer& _sc2;
};

#endif

