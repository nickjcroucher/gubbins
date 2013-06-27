// 	$Id: molphyFormat.cpp 962 2006-11-07 15:13:34Z privmane $	
#include "molphyFormat.h"
#include "someUtil.h"
#include "errorMsg.h"

sequenceContainer molphyFormat::read(istream &infile, const alphabet* alph) {
	sequenceContainer mySeqData = readUnAligned(infile, alph);
	mySeqData.makeSureAllSeqAreSameLengthAndGetLen();
	return mySeqData;
}
sequenceContainer molphyFormat::readUnAligned(istream &infile, const alphabet* alph) {

	vector<string> seqFileData;
	putFileIntoVectorStringArray(infile,seqFileData);
	if (seqFileData.empty()){
		errorMsg::reportError("unable to open file, or file is empty in molphy format");
	}
	vector<string>::iterator currentLinePosition = seqFileData.begin();

	string::const_iterator itStr = seqFileData.begin()->begin();
	string::const_iterator itStrEnd = seqFileData.begin()->end();

	int f_numSeq;
	bool readSeqNum= fromStringIterToInt(itStr,itStrEnd,f_numSeq);
	if (readSeqNum == false) errorMsg::reportError("Error reading number of sequences while reading MOLPHY sequence format");
	int f_seqLength;
	bool readSeqLen= fromStringIterToInt(itStr,itStrEnd,f_seqLength);
	if (readSeqLen == false) errorMsg::reportError("Error reading the sequences length while reading MOLPHY sequence format");
	currentLinePosition++; // we read the first line.

//---------------------------------------------------------------------
	sequenceContainer mySeqData;

//---------------------------------------------------------------------
//	vector<sequenceContainer::sequenceDatum*> vec;
//	seqDataPtr->getSequenceDatumPtrVectorNonConst(vec);

	int localID=-1;

	vector<string>::const_iterator it1 = seqFileData.begin();
	++it1; //skipping the first line that was read already.
	while (it1!= seqFileData.end()) {
		localID++;	  
		if (it1->empty()) {
			it1++;
			continue; // empty line continue
		}
		// read the name.
		string name(*it1);
		it1++;

		string tmpString;
		while (it1 != seqFileData.end()) {
			if (tmpString.size() < f_seqLength) {
				tmpString+=*it1;
				++it1;
			}
			else break;
		}
		
		mySeqData.add(sequence(tmpString,name,"",localID,alph));

	}
	return mySeqData;
}




void molphyFormat::write(ostream &out, const sequenceContainer& sd) {
	out<<sd.numberOfSeqs()<<" "<<sd.seqLen()<<endl;
	for (sequenceContainer::constTaxaIterator it5=sd.constTaxaBegin();it5!=sd.constTaxaEnd();++it5) {
		out<<it5->name()<<endl;
		string seqString = it5->toString();
		int k=0;
		for (string::const_iterator cPos=seqString.begin() ; cPos != seqString.end() ; cPos ++,k++ ) {
			if (k>0 && ((k%60)==0)) out<<endl;
			out<<*cPos;
		}
		out<<endl;
	}
}



