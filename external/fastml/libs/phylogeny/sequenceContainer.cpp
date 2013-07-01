// $Id: sequenceContainer.cpp 5244 2008-11-16 17:21:57Z cohenofi $
#include "sequenceContainer.h"
#include "logFile.h"
#include "someUtil.h"

sequenceContainer::sequenceContainer(const sequenceContainer& other,const alphabet *inAlph) :
_generalRemarks(other._generalRemarks),
_id2place(other._id2place)
{
	for (int i=0; i < other._seqDataVec.size(); ++i)
		_seqDataVec.push_back(sequence(other._seqDataVec[i],inAlph));
}


//if bAugumentShorterSeqs=true then add gap characters at the end of short seqeunces
const int sequenceContainer::makeSureAllSeqAreSameLengthAndGetLen(bool bAugumentShorterSeqs) {
	bAugumentShorterSeqs = true;
	if (_seqDataVec.size() == 0) return 0;
	const int len = _seqDataVec[0].seqLen();
	for (int i=1; i < _seqDataVec.size(); ++i) {
		if (_seqDataVec[i].seqLen()!=len) {
			if (bAugumentShorterSeqs) {
				for (int pos = _seqDataVec[i].seqLen(); pos < len; ++pos) 
					_seqDataVec[i].push_back(getAlphabet()->gap());
			}
			else {
                cerr<<_seqDataVec[i].name()<<" "<<_seqDataVec[i].seqLen()<<" "<<len<<endl;
                //errorMsg::reportError("not all sequences are of the same lengths");
			}
		}
	}

	return len;
}

//void sequenceContainer::addFromsequenceContainer(sequenceContainer& seqToAdd){
//	if (_seqDataVec.empty()) { // first sequence to add
//		sequenceContainer::taxaIterator tit;
//		sequenceContainer::taxaIterator titEND;
//		tit.begin(seqToAdd);
//		titEND.end(seqToAdd);
//		while (tit!=titEND) {
//			_seqDataVec.push_back(*tit);
//
//		}
//	}
//	else {// now we are adding sequences to sequences that are already there.
//		sequenceContainer::taxaIterator tit;
//		sequenceContainer::taxaIterator titEND;
//		tit.begin(seqToAdd);
//		titEND.end(seqToAdd);
//		while (tit!=titEND) {
//			for (int i=0; i < _seqDataVec.size(); ++i) {
//				if (tit->name() == _seqDataVec[i].name()) {
//					_seqDataVec[i]+=(*tit);
//					break;
//				}
//			}
//			++tit;
//		}
//	}
//}

void sequenceContainer::changeGaps2MissingData() {

	for (int i = 0; i < seqLen();++i) {//going over al positions
		for (int j = 0; j < _seqDataVec.size();++j) {
			if (_seqDataVec[j][i] == -1){
				 _seqDataVec[j][i]=getAlphabet()->unknown(); // missing data
			}
		}
	}
}

const int sequenceContainer::getId(const string &seqName, bool issueWarningIfNotFound) const {
	int k;
	for (k=0 ; k < _seqDataVec.size() ; ++k) {
		if (_seqDataVec[k].name() == seqName) return (_seqDataVec[k].id());
	}
	if (k == _seqDataVec.size() && issueWarningIfNotFound) {
		// debuggin
		LOG(5,<<"seqName = "<<seqName<<endl);
		for (k=0 ; k < _seqDataVec.size() ; ++k) {
			LOG(5,<<"_seqDataVec["<<k<<"].name() ="<<_seqDataVec[k].name()<<endl);
		}
		//end dubug
		LOG(0,<<seqName<<endl);
		vector<string> err;
		err.push_back("Could not find a sequence that matches the sequence name  ");
		err.push_back(seqName);
		err.push_back("in function sequenceContainer::getSeqPtr ");
		err.push_back(" make sure that names in tree file match name in sequence file ");
		errorMsg::reportError(err); // also quit the program
	}
	return -1;
}

const Vstring sequenceContainer::names() const {
	vector<string> res;
	for (int i=0; i < _seqDataVec.size(); ++i) {
		res.push_back(_seqDataVec[i].name());
	}
	return res;
}

sequenceContainer::sequenceContainer() {
	_id2place.resize(100,-1);
}

sequenceContainer::~sequenceContainer(){}

void sequenceContainer::add(const sequence& inSeq) {
	_seqDataVec.push_back(inSeq);
	if (_id2place.size() < inSeq.id()+1) {
		_id2place.resize(inSeq.id()+100,-1);
	}
	if (_id2place[inSeq.id()] != -1) {
		string err = "Two sequences with the same id - error in function sequenceContainer::add";
		err+= "\nThe id of the sequence you are trying to add = ";
		err += int2string(inSeq.id());
		errorMsg::reportError(err);
	}
	_id2place[inSeq.id()] = _seqDataVec.size()-1;
}


//given a sequence id the sequence is removed from the sequence container
//and the vector _id2place is updated.
void sequenceContainer::remove(const int idSeq)  {
	if (idSeq > _id2place.size()-1 || idSeq<0) 	
		errorMsg::reportError("the id of sequence is not mapped by id2place in function sequenceContainer::remove");
	int place = _id2place[idSeq];
	
	if (place < 0) 
		errorMsg::reportError("cannot find place of the id in the sequence container in function sequenceContainer::remove");
	_seqDataVec.erase(_seqDataVec.begin()+place);

	_id2place[idSeq] = -1;
	for (int i=place;i<_seqDataVec.size();i++) {
		int id = _seqDataVec[i].id();
		_id2place[id]--;
	}
}

 
//removes identical sequences in the sequence container.
void sequenceContainer::removeIdenticalSequences(){
	bool exist;
	for (int i=1;i<_seqDataVec.size();i++){
		sequence sq1 = _seqDataVec[i];
		for (int j=0;j<i;j++){
			sequence sq2 = _seqDataVec[j];
			exist = true;
			if (sq1.seqLen() != sq2.seqLen()) continue;
			for (int pos=0;pos<sq1.seqLen();pos++){
				if (sq1[pos] != sq2[pos]){
					exist = false;
					break;
				}
			}
			if (exist) { 
				remove(sq1.id());
				i--;
				break;
				
			}

		}
	
	}

}

void sequenceContainer::removeGapPositions(){
	vector<int> posToRemove(seqLen(),0);
	bool gapCol;
	int i,j;
	for (i = 0; i < seqLen();++i) {//going over al positions
		gapCol = false;
		for (j = 0; j < _seqDataVec.size();++j) {
			if (_seqDataVec[j][i] == -1) posToRemove[i] = 1;
		}
	}
	removePositions(posToRemove);
}
void sequenceContainer::removeGapPositionsAllSeqs(){
	vector<int> posToRemove(seqLen(),1);
	bool gapCol;
	int i,j;
	for (i = 0; i < seqLen();++i) {//going over al positions
		gapCol = false;
		for (j = 0; j < _seqDataVec.size();++j) {
			if (_seqDataVec[j][i] != -1) posToRemove[i] = 0;
		}
	}
	removePositions(posToRemove);
}
void sequenceContainer::removeGapPositionsAccordingToAReferenceSeq(const string & seqName){
	int idOfRefSeq = getId(seqName,true);
	vector<int> posToRemove(seqLen(),0);
	int i;
	for (i = 0; i < seqLen();++i) {//going over al positions
		if (_seqDataVec[idOfRefSeq][i] == -1) posToRemove[i] = 1;
	}
	removePositions(posToRemove);
}

void sequenceContainer::removeUnknownPositionsAccordingToAReferenceSeq(const string & seqName){
	int idOfRefSeq = getId(seqName,true);
	vector<int> posToRemove(seqLen(),0);
	int i;
	for (i = 0; i < seqLen();++i) {//going over al positions
		if (_seqDataVec[idOfRefSeq][i] == getAlphabet()->unknown()) posToRemove[i] = 1;
	}
	removePositions(posToRemove);
}

//removePositions: the positions to be removed are marked as '1' in posToRemoveVec
//all othehr positions are '0' 	
void sequenceContainer::removePositions(const Vint & posToRemoveVec) {
	for (int z = 0; z < _seqDataVec.size();++z) {
		_seqDataVec[z].removePositions(posToRemoveVec);
	}
}

void sequenceContainer::changeDotsToGoodCharacters() {
	for (int i = 0; i < seqLen();++i) {//going over al positions
		int charInFirstSeq = _seqDataVec[0][i];
		if (charInFirstSeq == -3) {
			LOG(5,<<" position is "<<i<<endl);
			errorMsg::reportError(" the first line contains dots ");
		}
		for (int j = 1; j < _seqDataVec.size();++j) {
			if ((_seqDataVec[j][i] == -3)) {
				_seqDataVec[j][i] = charInFirstSeq; // missing data
			}
		}
	}
}

int sequenceContainer::numberOfSequencesWithoutGaps (const int pos) const {
	int numOfNonCharPos = numberOfSeqs();
	for (int i=0; i < numberOfSeqs(); ++i) {
		if ((*this)[i][pos] <0) --numOfNonCharPos;
	}
	return numOfNonCharPos;
}

int sequenceContainer::numberOfSequencesWithoutUnknowns (const int pos) const {
	int numOfNonCharPos = numberOfSeqs();
	int unknown = getAlphabet()->unknown();
	for (int i=0; i < numberOfSeqs(); ++i) {
		if ((*this)[i][pos] == unknown ) 
			--numOfNonCharPos;
	}
	return numOfNonCharPos;
}

bool sequenceContainer::isInvariable(const int pos) const {
	int charFound = getAlphabet()->unknown(); 
	for (int i=0; i < numberOfSeqs(); ++i) {
		if ((*this)[i][pos] >= 0) {
			if (charFound == getAlphabet()->unknown())
				charFound = (*this)[i][pos];
			else if (charFound != (*this)[i][pos])
				return false;
		}
	}
	return true;
}

int sequenceContainer::getInvariablePosNum() const {
	int sum = 0;
	for (int pos = 0; pos < seqLen(); ++pos) {
		if (isInvariable(pos))
			++sum;
	}
	return sum;
}

// new func for gainLoss project
void sequenceContainer::startZeroSequenceContainerGL(const sequenceContainer &sc, const gainLossAlphabet& alph, const int minNumOfOnes)
{
	//if(minNumOfOnes == 1){
	//	string str = "0";
	//	string remark;
	//	int localid =0;	
	//	for(int i=0; i<sc.numberOfSeqs();i++){
	//		this->add(sequence(str,sc.name(i),remark,localid,&alph));
	//		++localid;
	//	}
	//}
	string str0 = "0";
	string str1 = "1";
	vector<string> strV;
	strV.resize(sc.numberOfSeqs());
	string remark ="";
	switch (minNumOfOnes) {
			case (1) :
				for(int i=0; i<sc.numberOfSeqs();i++){
					// add patterns of 0 ones
					strV[i] = str0;
				}
				break;
			case (2) :
				for(int i=0; i<sc.numberOfSeqs();i++){
					// add patterns of 0 ones
					strV[i] = str0;
				}
				for(int i=0; i<sc.numberOfSeqs();i++){
					// add patterns of only 1 ones
					for(int j=0; j<sc.numberOfSeqs(); j++){
						if(j==i){
							strV[i]+=str1;
						}
						else{
							strV[i]+=str0;
						}
					}
				}
				break;
			case (3) :
				for(int i=0; i<sc.numberOfSeqs();i++){
					// add patterns of 0 ones
					strV[i] = str0;
				}
				for(int i=0; i<sc.numberOfSeqs();i++){
					// add patterns of only 1 ones
					for(int j=0; j<sc.numberOfSeqs(); j++){
						if(j==i){
							strV[i]+=str1;
						}
						else{
							strV[i]+=str0;
						}
					}
				}
				// add patterns of only 2 ones
				for(int onePosition1=0; onePosition1<sc.numberOfSeqs(); onePosition1++){
					for(int onePosition2=0; onePosition2<sc.numberOfSeqs(); onePosition2++){
						if(onePosition2<=onePosition1)
							continue;
						for(int i=0; i<sc.numberOfSeqs();i++){
							if(i==onePosition1 || i==onePosition2){
								strV[i]+=str1;
							}
							else{
								strV[i]+=str0;
							}			
						}
					}
				}				
				break;
	}
	for(int i=0; i<sc.numberOfSeqs();i++){
		//cout<<strV[i]<<endl;
		this->add(sequence(strV[i],sc.name(i),remark,i,&alph));
	}
}



//concatenate two sequecneContainers. 
//The sequence names must be identical in the two containers.
//returns false if: (1) A sequence_name in one of the containers does not match any sequence_name in the other container.
bool sequenceContainer::concatenate(const sequenceContainer& other) {
	if (other.numberOfSeqs() != numberOfSeqs())
		return false;
	for(int i = 0; i < numberOfSeqs(); ++i)
	{
		bool bFound = false;
        for (int j = 0; j < other.numberOfSeqs(); ++j)
		{
		    if((*this)[i].name() == other[j].name())
			{
				(*this)[i] += other[i];
				bFound = true;
				break;
			}
		}
		if (bFound == false) 
		{
			string msg = string("Can't find sequence name in the second MSA: ") + other[i].name();
            errorMsg::reportError(msg);
			return false;
		}
	}
	return true;
}
