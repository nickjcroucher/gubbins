// $Id: amino.h 1901 2007-03-15 13:21:06Z nimrodru $

#ifndef ____AMINO
#define ____AMINO

#include "definitions.h"
#include "errorMsg.h"
#include "alphabet.h"
#include "geneticCodeHolder.h"
#include "codon.h"


//utility of amino acid
class aminoUtility {
public:

	static vector<int> codonOf(const int a, codon &cod); //returns vector of codons that code to a under a specific genetic code.
	
};

//based on the amino-acid list found in http://www.dur.ac.uk/~dbl0www/Bioinformatics/aminoacids.htm
class amino : public alphabet {
public:
	explicit amino();
	virtual ~amino() {}
	virtual alphabet* clone() const { return new amino(*this); }
	int unknown() const  {return -2;}
	int gap() const  {return -1;}
	int size() const {return 20;}
	int stringSize() const {return 1;} // one letter code.
	int relations(const int charInSeq, const int charToCheck) const;
	int fromChar(const string& str, const int pos) const;
	int fromChar(const char s) const;
	string fromInt(const int in_id) const;
	vector<int> fromString(const string& str) const;
	// "specific" here is not unknown, nor ambiguity, nor gap (for example, for nucleotides it will true for A,C,G, or T).
	bool isSpecific(const int id) const {return (id>=0 && id < size());}

private:
	int relations_internal(const int charInSeq, const int charToCheck) const;
	VVint _relation;
};//end of class

#endif


