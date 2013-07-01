// $Id: sequenceContainer.h 5244 2008-11-16 17:21:57Z cohenofi $

#ifndef ___SEQUENCE_CONTAINER
#define ___SEQUENCE_CONTAINER
#include "definitions.h"
#include "sequence.h"
#include "gainLossAlphabet.h"

class sequenceContainer {
public:

	class taxaIterator;
	friend class taxaIterator;
	class constTaxaIterator;
	friend class constTaxaIterator;

//------------------------------------------------------------
//constructors:
    explicit sequenceContainer();
	sequenceContainer(const sequenceContainer& other,const alphabet *inAlph);
	virtual ~sequenceContainer();

	//questions only:
	const int seqLen() const {return _seqDataVec.empty()? 0 : _seqDataVec[0].seqLen();}
	const int numberOfSeqs() const {return _seqDataVec.size();}
	const int alphabetSize() const {return _seqDataVec.empty()? 0 : _seqDataVec[0].getAlphabet()->size();}
	const vector<string>& getGeneralRemarks() const {return _generalRemarks;}
	const int makeSureAllSeqAreSameLengthAndGetLen(bool bAugumentShorterSeqs = false); //if bAugumentShorterSeqs=true then add gap characters at the end of short seqeunces
	const int getId(const string &seqName, bool issueWarninInNotFound=true) const;//return -1 if not found...
	sequence& operator[](const int id) {return  _seqDataVec[_id2place[id]];} // get the ID of the sequence. Return the sequence itself.
	const sequence& operator[](const int id) const {return _seqDataVec[_id2place[id]];}
	const Vstring names() const; // return a vector<string> of the names of all the sequences.
	const string& name(const int id) const {return _seqDataVec[_id2place[id]].name();};
	const alphabet* getAlphabet() const {return _seqDataVec[0].getAlphabet();}
	//returns the number of positions that are invariable (all seqs are identical
	int getInvariablePosNum() const;
	bool isInvariable(const int pos) const; 
	// computed the number of sequences without gaps at a specific position
	// for example, if the multiple sequence alignment is 
	// AT-
	// AG-
	// A-M
	// numberOfSequencesWithoutGaps(0) = 3
	// numberOfSequencesWithoutGaps(1) = 2
	// numberOfSequencesWithoutGaps(2) = 1
	int numberOfSequencesWithoutGaps(const int pos) const;
	int numberOfSequencesWithoutUnknowns(const int pos) const;
	
	
//make changes:
	void resize(int t,const alphabet* inAlph) {
		if (inAlph == NULL) {
			errorMsg::reportError("cannot resize when the alphabet is unknown");
		}
		sequence s(inAlph);
		_seqDataVec.resize(t,s);
	}
	void add(const sequence& inSeq);
	void remove(const int idSeq);
	void removeIdenticalSequences();
	int placeToId(const int place) const {return _seqDataVec[place].id();}; //get place in the vector and return the id of the sequence
	void addGeneralRemark(const string& inRemark) {_generalRemarks.push_back(inRemark);}
	void changeGaps2MissingData();
	//removePositions: the positions to be removed are marked as '1' in posToRemoveVec
	//all othehr positions are '0' 	
	void removePositions(const Vint & posToRemoveVec);	
	void removeGapPositions();
	void removeGapPositionsAllSeqs();
	void removeGapPositionsAccordingToAReferenceSeq(const string & seqName);
	void changeDotsToGoodCharacters();
	void removeUnknownPositionsAccordingToAReferenceSeq(const string & seqName);
	bool concatenate(const sequenceContainer& other);
	void startZeroSequenceContainerGL(const sequenceContainer &sc, const gainLossAlphabet& alph, const int minNumOfOnes=1);
	

public: 
	sequence::Iterator begin(const int id){//iterface to sequence iterator
		sequence::Iterator temp;
		temp.begin(_seqDataVec[id]);
		return temp;
	}
	sequence::Iterator end(const int id){//iterface to sequence iterator
		sequence::Iterator temp;
		temp.end(_seqDataVec[id]);
		return temp;
	}

	class taxaIterator {
	public:
		explicit taxaIterator(){};
		~taxaIterator(){};
		void begin(sequenceContainer & inSeqCont){
			_pointer = inSeqCont._seqDataVec.begin();
		}
	    void end(sequenceContainer & inSeqCont){
			_pointer = inSeqCont._seqDataVec.end();
		}
		sequence& operator* ()  {return *_pointer;}
		sequence const &  operator* () const {return *_pointer;}
		sequence *  operator-> ()  {return &*_pointer;} //MATAN- CHECK!!!
		sequence const *  operator-> () const {return &* _pointer;} // MATAN - CHECK!!!

		void operator ++() {++_pointer;}
	    void operator --() { --_pointer; }
	    bool operator != (const taxaIterator& rhs){return (_pointer != rhs._pointer);}
	    bool operator == (const taxaIterator& rhs){return (_pointer == rhs._pointer);}
	private:
		vector<sequence>::iterator _pointer;
	};//end if class taxaIterator


	class constTaxaIterator {
	public:
		explicit constTaxaIterator(){};
		~constTaxaIterator(){};
	    void begin(const sequenceContainer & inSeqCont){
			_pointer = inSeqCont._seqDataVec.begin();
		}
		void end(const sequenceContainer & inSeqCont){
			_pointer = inSeqCont._seqDataVec.end();
		}
		sequence const &  operator*() const {return *_pointer;}
		sequence const *  operator->() const {return &*_pointer;}// MATAN - CHECK!!!

		void operator ++() {++_pointer;}
		void operator --() { --_pointer; }
		bool operator != (const constTaxaIterator& rhs) {
		  return (_pointer != rhs._pointer);
		}

		bool operator == (const constTaxaIterator& rhs) {
		  return (_pointer == rhs._pointer);
		}
	private:
		vector<sequence>::const_iterator _pointer;
	};

	public: // interfaces to iterators
	taxaIterator taxaBegin(const int id=0){// interface to taxaIterator
		taxaIterator temp;
		temp.begin(*this);
		return temp;
	}

	taxaIterator taxaEnd(){// interface to taxaIterator
		taxaIterator temp;
		temp.end(*this);
		return temp;
	}

	constTaxaIterator constTaxaBegin() const{ //interface to const taxaIter
		constTaxaIterator temp;
		temp.begin(*this);
		return temp;
	}
	constTaxaIterator constTaxaEnd() const{
		constTaxaIterator temp;
		temp.end(*this);
		return temp;
	  }

	private:
	vector<sequence> _seqDataVec;
	vector<string> _generalRemarks;
	vector<int> _id2place;
};

#endif

