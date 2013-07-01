// $Id: fastaFormat.cpp 962 2006-11-07 15:13:34Z privmane $
#include "fastaFormat.h"
#include "someUtil.h"
#include "errorMsg.h"
#include <algorithm>
using namespace std;

sequenceContainer fastaFormat::read(istream &infile, const alphabet* alph) {
	sequenceContainer mySeqData = readUnAligned(infile, alph);
	mySeqData.makeSureAllSeqAreSameLengthAndGetLen();
	return mySeqData;
}


sequenceContainer fastaFormat::readUnAligned(istream &infile, const alphabet* alph) {
	sequenceContainer mySeqData;

	vector<string> seqFileData;
	putFileIntoVectorStringArray(infile,seqFileData);
	if (seqFileData.empty()){
		errorMsg::reportError("unable to open file, or file is empty in fasta format");
	}

	vector<string>::const_iterator it1;
	int localid=0;
	for (it1 = seqFileData.begin(); it1!= seqFileData.end(); ) {
		if (it1->empty()) {++it1;continue; }// empty line continue

		string remark;
		string name;

		if ((*it1)[0] == '>') {
			string::const_iterator itstrtmp = (*it1).begin();
			itstrtmp++;
			while (itstrtmp != (*it1).end()) {
				name+= *itstrtmp;
				itstrtmp++;
			}

			//for (string::iterator i = name.begin(); i!=(name.end()-2);++i) {
			//	*i=*(i+1); // removing the ">". should be done more elegant...
			//}
			++it1;
		} else {
			LOG(0,<<"problem in line: "<<*it1<<endl);
			errorMsg::reportError("Error reading fasta file, error finding sequence name starting with >",1);
		}
		while (it1->empty()) it1++; // empty line continue
		
		string str;
		while (it1!= seqFileData.end()) {
			if ((*it1)[0] == '>') break;
			str+=*it1;
			++it1;
		}
		// remove spaces form str;
		str.erase(
			std::remove(str.begin(),str.end(),' '),str.end()
			);

		mySeqData.add(sequence(str,name,remark,localid,alph));
		localid++;
	}

	return mySeqData;
}


void fastaFormat::write(ostream &out, const sequenceContainer& sd) {
	for (sequenceContainer::constTaxaIterator it5=sd.constTaxaBegin();it5!=sd.constTaxaEnd();++it5) {
		out<<">"<<(it5)->name()<<endl;
		out<<it5->toString()<<endl;
	}
}

