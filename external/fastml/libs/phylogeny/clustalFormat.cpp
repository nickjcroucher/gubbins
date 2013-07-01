// $Id: clustalFormat.cpp 962 2006-11-07 15:13:34Z privmane $

#include "clustalFormat.h"
#include "codon.h"
#include "someUtil.h"
#include "errorMsg.h"
#include <map>

sequenceContainer clustalFormat::read(istream &infile, const alphabet* alph) {
	sequenceContainer mySeqData = readUnAligned(infile, alph);
	mySeqData.makeSureAllSeqAreSameLengthAndGetLen();
	return mySeqData;
}

sequenceContainer clustalFormat::readUnAligned(istream &infile, const alphabet* alph) {
	sequenceContainer mySequenceData;

	vector<string> seqFileData;
	map<string ,string> stringsToAdd;   //map that holding for each name last
										//one or two nucleotides (when reading codon
										//alphabet) of the line in order to add it
										//to the next line.
	putFileIntoVectorStringArray(infile,seqFileData);
	if (seqFileData.empty()){
		errorMsg::reportError("unable to open file, or file is empty in clustal format");
	}


	vector<string>::const_iterator it1= seqFileData.begin();

	// make sure that the first 7 chars in the first line is clustal
	if (it1->size()<7) errorMsg::reportError("first word in clusltal sequence file format must be clustal",1);
	if (  (( (*it1)[0] != 'C') && ((*it1)[0] != 'c'))
		|| (((*it1)[1] != 'L') && ((*it1)[1] != 'l'))
		|| (((*it1)[2] != 'U') && ((*it1)[2] != 'u'))
		|| (((*it1)[3] != 'S') && ((*it1)[3] != 's'))
		|| (((*it1)[4] != 'T') && ((*it1)[4] != 't'))
		|| (((*it1)[5] != 'A') && ((*it1)[5] != 'a'))
		|| (((*it1)[6] != 'L') && ((*it1)[6] != 'l')) ) {
		errorMsg::reportError("first word in clusltal sequence file format must be clustal",1);
	}
	it1++;

	int localid=0;
	while (it1!= seqFileData.end()) {
		if (it1->empty()) {++it1;continue; }// empty line continue
		if ((it1->size() > 1) && ((*it1)[0]==' ')) {++it1;continue; }// remark line 
		string remark;
		string name;

//		getFromLineAnameAndAsequence;
		string name1;
		string stringSeq1;
		string::const_iterator it2 = (it1)->begin();
		for (; it2 != (it1)->end();++it2) {
			if ((*it2)==' ') break;
			else name1+=(*it2);
		}
		if (stringsToAdd.find(name1)!=stringsToAdd.end()) //not new sequence
			stringSeq1 = stringsToAdd[name1]; //init stringSeq1 with the nucleotide 
											  //from the previous line
		for (; it2 != (it1)->end();++it2) {
			if ((*it2)==' ') continue;
			else stringSeq1+=(*it2);
		}
		
		//when alphabet is codon stringSeq1 must be product of three.
		// 1. save  1 or  2 last nucleotide in stringToAdd
		// 2. substr the last or two last nucleotide for the next line.
		// 3. keep stringToAdd in map (according the name).
		string stringToAdd="";
		//		codon codonAlph;
		if (alph->size()>=60){	// codon?
			if ((stringSeq1.size()%3)==1){  //add the last nucleotide to the next line
				stringToAdd+=stringSeq1[stringSeq1.size()-1];
				stringSeq1 = stringSeq1.substr(0,stringSeq1.size()-1);
			}
			if ((stringSeq1.size()%3)==2){ //add the 2 last nucleotide to the next line
				stringToAdd+=stringSeq1[stringSeq1.size()-2];
				stringToAdd+=stringSeq1[stringSeq1.size()-1];
				stringSeq1 = stringSeq1.substr(0,stringSeq1.size()-2);
			}
		
		}
		stringsToAdd[name1] = stringToAdd; //update the map with the new stringToAdd 
		int id = mySequenceData.getId(name1,false);
		if (id==-1) { // new sequence.
			name = name1;
			mySequenceData.add(sequence(stringSeq1,name,remark,localid,alph));
			localid++;
		} else {// the sequence is already there...
			sequence tmp(stringSeq1,name,remark,id,alph);
			mySequenceData[id].operator += (tmp);
		}
	
		it1++;
	}

	return mySequenceData;
}

void clustalFormat::write(ostream &out, const sequenceContainer& sd) {
	// setting some parameters
	const int numOfPositionInLine = 60;
	int maxLengthOfSeqName =0;
	for (sequenceContainer::constTaxaIterator p=sd.constTaxaBegin(); p != sd.constTaxaEnd(); ++p ) {
		int nameLen = (*p).name().size();
		if (nameLen>maxLengthOfSeqName) maxLengthOfSeqName=nameLen;
	}
	if (maxLengthOfSeqName<15) maxLengthOfSeqName=16;
	else maxLengthOfSeqName=maxLengthOfSeqName+4; // all this maxLengthOfSeqName is the 

	out<<"CLUSTAL V"<<endl;
	           // num. of space after the name.
	int currentPosition = 0;
	int charLen = sd.seqLen();
	//in case of codon alphabet the character length is : 3*(sequence_length)
	//	codon codonAlph;
	if (sd.alphabetSize()>=60) charLen*=3; 
	out<<endl<<endl;
	while (currentPosition < charLen ) {
		out.flush();
		//for (vector<const sequenceContainer::sequenceDatum*>::const_iterator it5= vec.begin(); it5!=vec.end(); ++ it5) {
		for (sequenceContainer::constTaxaIterator it5=sd.constTaxaBegin();it5!=sd.constTaxaEnd();++it5) {
			for (int iName = 0 ;iName<maxLengthOfSeqName; ++iName) {
				if (iName<(*it5).name().size()) {
					out<<(*it5).name()[iName];
					out.flush();
				}
				else out<<" ";
			}
			out.flush();
			out<<" ";
			
			if (charLen<numOfPositionInLine) 
				out<<it5->toString()<<endl;
			else {
				for (int k=currentPosition; k < currentPosition+numOfPositionInLine; ++k) {
					if (k>=charLen) 
						break;
					out<<it5->toString()[k];
					//in case of codon alphabet each position is three characters
					
					if (sd.alphabetSize()>=60){ 
						out<<it5->toString()[++k];
						out<<it5->toString()[++k];
					}
				}
				out<<endl;
			}
		}
		currentPosition +=numOfPositionInLine;
		out<<endl<<endl;
	}

	return;
}

