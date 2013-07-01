// $Id: phylipFormat.cpp 962 2006-11-07 15:13:34Z privmane $

#include "phylipFormat.h"
#include "someUtil.h"
#include "errorMsg.h"
#include "logFile.h"

sequenceContainer phylipFormat::read(istream &infile, const alphabet* alph){
	sequenceContainer mySeqData = readUnAligned(infile, alph);
	mySeqData.makeSureAllSeqAreSameLengthAndGetLen();
	return mySeqData;
}
sequenceContainer phylipFormat::readUnAligned(istream &infile, const alphabet* alph){
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
		if (currentLinePosition->empty()) {++currentLinePosition;continue;} // empty line constinue
		string remark;
		string name;
		sequence seq(alph);


		
		if (mySeqData.numberOfSeqs() < f_numSeq ) {//get from the line a name and a sequence;
			
			string name1;
			string stringSeq1;
			string::const_iterator it2 = (currentLinePosition)->begin();
			for (; it2 != (currentLinePosition)->end();++it2) {
				if ((*it2)==' ') break;
				else name1+=(*it2);
			}
			for (; it2 != (currentLinePosition)->end();++it2) {
				if ((*it2)==' ') continue;
				else stringSeq1+=(*it2);
			}
			mySeqData.add(sequence(stringSeq1,name1,remark,localid,alph));
			currentLinePosition++;
			localid++;
		}
		else { // adding to the 
			string stringSeq1;
			string::const_iterator it2 = (currentLinePosition)->begin();
			int sequenceId=localid%f_numSeq;
			for (; it2 != (currentLinePosition)->end() && 
			       mySeqData[sequenceId].seqLen() <f_seqLength;++it2) {
				if ((*it2)==' ') continue;
				else stringSeq1+=(*it2);
				
			}
			sequence tmp(stringSeq1,"","",sequenceId,alph);
			mySeqData[sequenceId].operator += (tmp);
			currentLinePosition++;
			localid++;
		}
	}
	return mySeqData;
}

void phylipFormat::write(ostream &out, const sequenceContainer& sd,
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

	int currentPosition = 0;
	while (currentPosition < sd.seqLen() ) {
		out<<endl;
		out.flush();
		// for (vector<const sequenceContainer::sequenceDatum*>::const_iterator it5= vec.begin(); it5!=vec.end(); ++ it5) {
		   for (sequenceContainer::constTaxaIterator it5=sd.constTaxaBegin();it5!=sd.constTaxaEnd();++it5) {

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
			
			if (it5->seqLen()<numOfPositionInLine) 
				out<<it5->toString()<<endl;
			else {
				for (int k=currentPosition; k < currentPosition+numOfPositionInLine; ++k) {
					if (k>=it5->seqLen()) break;
					out<<it5->toString(k);
					if (((k+1)%spaceEvery==0) && (((k+1)%numOfPositionInLine!=0))) out<<" ";
				}
				out<<endl;
			}
		}
		currentPosition +=numOfPositionInLine;
		
	}
	return;
}


