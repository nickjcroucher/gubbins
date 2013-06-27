// $Id: mulAlphabet.h 1901 2007-03-15 13:21:06Z nimrodru $

// version 1.01
// last modified 1 Jan 2004

#ifndef ___MUL_ALPHABET_H
#define ___MUL_ALPHABET_H

#include "definitions.h"
#include "alphabet.h"
#include "someUtil.h"

class mulAlphabet : public alphabet {

public:
	mulAlphabet(const alphabet* baseAlphabet, int mulFactor);
	mulAlphabet(const mulAlphabet& other);
	virtual ~mulAlphabet();
    virtual alphabet* clone() const { return new mulAlphabet(*this); }
	mulAlphabet& operator=(const mulAlphabet &other);

	int unknown() const ;
	int gap() const;
	
	int size() const {return _size;}
	int stringSize() const ; 
	bool isSpecific(const int id) const ;

	int fromChar(const string& str, const int pos) const;
	vector<int> fromString(const string& str) const;

	string fromInt(const int id) const;

	int relations(const int charInSeq, const int charToCheck) const;
	int compareCategories(int charA, int charB) const;
	const alphabet* getBaseAlphabet() const {return _baseAlphabet;}
	
public:
	int convertFromBasedAlphaInt(int id) const;
	int convertToBasedAlphaInt(int id) const;	
	
private:
	alphabet* _baseAlphabet; // This alphabet must use single characters, i.e. - not codon. (or we will have to add to every alphabet a member which holds its character's size)
	int _mulFactor ; // number of times that the alphabet is multiplied by = Number of categories (g in Galtier paper)
	int _size ; // this is simply the _baseAlphabet->size() * _mulFactor


};

#endif

