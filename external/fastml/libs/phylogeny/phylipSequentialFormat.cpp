// $Id: phylipFormat.cpp 962 2006-11-07 15:13:34Z privmane $

#include "phylipSequentialFormat.h"
#include "someUtil.h"
#include "errorMsg.h"
#include "logFile.h"

sequenceContainer phylipSequentialFormat::read(istream &infile, const alphabet* alph){
	sequenceContainer mySeqData = readUnAligned(infile, alph);
	mySeqData.makeSureAllSeqAreSameLengthAndGetLen();
	return mySeqData;
}
sequenceContainer phylipSequentialFormat::readUnAligned(istream &infile, const alphabet* alph){
	sequenceContainer mySeqData;

	vector<string> seqFileData;
	putFileIntoVectorStringArray(infile,seqFileData);

	vector<string>::const_iterator currentLinePosition = seqFileData.begin();
	string::const_iterator itStr = seqFileData.begin()->begin();
	string::const_iterator itStrEnd = seqFileData.begin()->end();

	int f_numSeq;
	bool readSeqNum= fromStringIterToInt(itStr,itStrEnd,f_numSeq);
	if (readSeqNum == false) errorMsg::reportError("Error reading number of sequences while reading PHYLIP sequence format");
	int f_seqLength;
	bool readSeqLen= fromStringIterToInt(itStr,itStrEnd,f_seqLength);
	if (readSeqLen == false) errorMsg::reportError("Error reading the sequences length while reading PHYLIP sequence format");
	currentLinePosition++; // we read the first line.

	int localid=0;
	for (; currentLinePosition != seqFileData.end() ; ) {
		if (currentLinePosition->empty()) {++currentLinePosition;continue;} // empty line continue
		string stringSeq1;
		string name1;
		while (stringSeq1.length() < f_seqLength ) { // adding a new seq			
			string::const_iterator it2 = (currentLinePosition)->begin();
			if ((*it2)==' ') { // line without seq. name, read seq. content only
				for (; it2 != (currentLinePosition)->end();++it2) {
					if ((*it2)==' ') continue;
					else stringSeq1+=(*it2);
				}
			}
			else { // first read sequence name, then read seq itself
				for (; it2 != (currentLinePosition)->end();++it2) {
					if ((*it2)==' ') break;
					else name1+=(*it2);
				}
				for (; it2 != (currentLinePosition)->end();++it2) {
					if ((*it2)==' ') continue;
					else stringSeq1+=(*it2);
				}
			}
			
			currentLinePosition++;
		}
		mySeqData.add(sequence(stringSeq1,name1,"",localid,alph));
		localid++;

	}
	return mySeqData;
}

void phylipSequentialFormat::write(ostream &out, const sequenceContainer& sd,
						 const int numOfPositionInLine,
						 const int spaceEvery) {
	sequenceContainer::constTaxaIterator it5=sd.constTaxaBegin();
	for (;it5!=sd.constTaxaEnd();++it5) {
		if (it5->name().size() > 10) break;
	}
	if (it5 != sd.constTaxaEnd()) {
		LOG(1,<<"you asked to print in phylip format\n");
		LOG(1,<<"however, the names in phylip format\n");
		LOG(1,<<"must be no more than 10 characters.\n");
		LOG(1,<<"Names are hence trancated to ten   \n");
		LOG(1,<<"characters. Notice, that this might\n");
		LOG(1,<<"result in a two or more sequences  \n");
		LOG(1,<<"having the same name               \n");
	}
	
	//	vector<const sequenceContainer::sequenceDatum*> vec;
	//	sd.getSequenceDatumPtrVector(vec);
	out<<sd.numberOfSeqs()<<"   "<<sd.seqLen();
	if (sd.constTaxaBegin()==sd.constTaxaEnd()) return;
	
	int maxLengthOfSeqName =0;
	maxLengthOfSeqName=10;	// all this maxLengthOfSeqName is the 


	for (sequenceContainer::constTaxaIterator it5=sd.constTaxaBegin();it5!=sd.constTaxaEnd();++it5) {
		int currentPosition = 0;
		out<<endl;
		out.flush();
		// first - print name of sequence
		for (int iName = 0 ;iName<maxLengthOfSeqName; ++iName) {
			if (iName<it5->name().size()) {
				if (currentPosition<numOfPositionInLine) {
					out<<it5->name()[iName];
				}
				else out<<" ";
				out.flush();
			}
			else out<<" ";
		}
		out.flush();
		out<<" ";
		// next - print sequence itself
		while (currentPosition < sd.seqLen() ) {
			if (it5->seqLen()<numOfPositionInLine) 
				out<<it5->toString()<<endl;
			else {
				for (int k=currentPosition; k < currentPosition+numOfPositionInLine; ++k) {
					if (k>=it5->seqLen()) break;
					out<<it5->toString(k);
					if (((k+1)%spaceEvery==0) && (((k+1)%numOfPositionInLine!=0))) out<<" ";
				}
				out<<endl;
				if (currentPosition+numOfPositionInLine < sd.seqLen()) {
				for (int i = 0; i <  spaceEvery +1; i++) // creates spaces to align properly
					out << " ";
				}
			}
			currentPosition +=numOfPositionInLine;
		}
		
	}

}


