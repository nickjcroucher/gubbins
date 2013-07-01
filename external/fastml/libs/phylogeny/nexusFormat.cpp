// $Id: nexusFormat.cpp 5987 2009-03-18 18:13:53Z itaymay $

#include "nexusFormat.h"
#include "someUtil.h"
#include "errorMsg.h"
#include <map>

sequenceContainer nexusFormat::read(istream &infile, const alphabet* pAlph) {
	sequenceContainer mySeqData = readUnAligned(infile, pAlph);
	mySeqData.makeSureAllSeqAreSameLengthAndGetLen();
	return mySeqData;
}

sequenceContainer nexusFormat::readUnAligned(istream &infile, const alphabet* pAlph) {
	if (!infile) {
		errorMsg::reportError("unable to read mase format, could not open file");
	}
	sequenceContainer mySeqData;;

	vector<string> seqFileData;
	putFileIntoVectorStringArray(infile,seqFileData);

	vector<string>::const_iterator it1 = seqFileData.begin();
	// make sure that the first 6 chars in the first line is #NEXUS
	if (it1->size()<6) errorMsg::reportError("first word in a nexus sequence file format must be #NEXUS",1);
	if (  ((*it1)[0] != '#') 
		|| (((*it1)[1] != 'N') && ((*it1)[1] != 'n'))
		|| (((*it1)[2] != 'E') && ((*it1)[2] != 'e'))
		|| (((*it1)[3] != 'X') && ((*it1)[3] != 'x'))
		|| (((*it1)[4] != 'U') && ((*it1)[4] != 'u'))
		|| (((*it1)[5] != 'S') && ((*it1)[5] != 's')) ) {
		errorMsg::reportError("first word in a nexus sequence file format must be #NEXUS",1);
	}
	it1++;

	while ( ( (*it1).find("matrix")  == -1) && ( (*it1).find("MATRIX")  == -1) && (it1!= seqFileData.end()))
	{ //check for the word matrix
		++it1;
	}
	
	int localid=0;
	//int x1 = ((*it1).find("matrix") != -1);
	//int x2 = ((*it1).find("MATRIX") != -1);
	if (((*it1).find("matrix") != -1) || ((*it1).find("MATRIX") != -1))
	{
		//taken from clustalFormat:
        //In case of codon alpahabet we cannot add a seqeunce that is not dividable by 3.
		//In this case the last nucleotides in each line (zero, one or two) 
		//should be saved. The next time the same sequence name appears - 
		//these saveed nucleotidea and are added to the begining of the line.
		map<string ,string> stringsToAdd;   


		for (++it1; it1 != seqFileData.end() ; ++it1)
		{
			if (((*it1).find("end;") != -1) || ((*it1).find("END;") != -1))
				break;
			if (it1->empty() || ((*it1).find(';') != -1)) 
			{ // empty line constinue
				continue;
			}
			sequence seq(pAlph);
			
			string taxonName;
			string remark;
			string stringSeq;
			bool beforeName = true;
			string::const_iterator stringIt = (it1)->begin();
			for (; stringIt != (it1)->end(); ++stringIt)
			{ //first loop finds the taxon name
				if ( ((*stringIt) == ' ') || ((*stringIt) == '\t'))
					if (beforeName == true)
						continue; //spaces before taxon name are legal
					else
						break; //A space marks the end of the taxon name 
				else 
				{
					taxonName += (*stringIt);
					beforeName = false;
				}
			}

			//check if a new sequence.
			//if the name already exists then init stringSeq with the nucleotide from the previous line of the same sequence
			if (stringsToAdd.find(taxonName)!=stringsToAdd.end()) 
				stringSeq = stringsToAdd[taxonName]; 

			for (; stringIt != (it1)->end(); ++stringIt) 
			{//second loop finds the sequecne
				if ( ((*stringIt)==' ') || 	((*stringIt) == '\t'))
					continue;
				else stringSeq += (*stringIt);
			}
			
			//when alphabet is codon stringSeq must be dividable by 3.
			// 1. save the reminder (0,1 or 2 last nucleotides) in stringToAdd
			// 2. substr the reminder from the sequence line.
			// 3. keep stringToAdd in map (according the name) to be added later.
			string stringToAdd="";
			if (pAlph->size()>=60){	// codon?
				if ((stringSeq.size()%3)==1){  //add the last nucleotide to the next line
					stringToAdd += stringSeq[stringSeq.size()-1];
					stringSeq = stringSeq.substr(0,stringSeq.size()-1);
				}
				if ((stringSeq.size() % 3) == 2){ //add the 2 last nucleotide to the next line
					stringToAdd+=stringSeq[stringSeq.size()-2];
					stringToAdd+=stringSeq[stringSeq.size()-1];
					stringSeq = stringSeq.substr(0, stringSeq.size() - 2);
				}
			}
			stringsToAdd[taxonName] = stringToAdd; //update the map with the new stringToAdd 
			//add sequence to container
			int id = mySeqData.getId(taxonName, false);
            if (id==-1) { // new sequence.
                mySeqData.add(sequence(stringSeq, taxonName,remark,localid, pAlph));
                localid++;
			}
			else {// the sequence is already there...
				sequence tmp(stringSeq,taxonName, remark, id, pAlph);
				mySeqData[id].operator += (tmp);
			}
		}
	}
	else
	{
		errorMsg::reportError("no sequence data in nexus file - no matrix keyword found");
	}
	
	return mySeqData;
}

void nexusFormat::write(ostream &out, const sequenceContainer& sc) {
	//vector<string> gfr = sd.getGeneralFileRemarks();
	//if (gfr.empty()) out<<";;\n;;\n";
	//for (vector<string>::const_iterator k=gfr.begin() ; k !=  gfr.end() ; ++k )
	//	out<<(*k)<<endl;
	out<<"#NEXUS"<<endl;
	out<<"begin data;"<<endl;
	out<<"dimensions ntax="<<sc.numberOfSeqs()<<" nchar="<<sc.seqLen() <<";"<<endl;
	if (sc.alphabetSize() == 4) 
		out<<"format datatype=dna gap=-;"<<endl;
	else
		out<<"format datatype=protein gap=-;"<<endl;
	out<<"matrix"<<endl;

	for (sequenceContainer::constTaxaIterator itSeq=sc.constTaxaBegin();itSeq!=sc.constTaxaEnd();++itSeq) {
		out<<"\t"<<itSeq->name()<<"\t"<<itSeq->toString()<<endl;
	}
	out<<";"<<endl;
	out<<"end;"<<endl;
}

