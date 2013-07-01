// $Id: mulAlphabet.cpp 1927 2007-04-04 16:44:23Z privmane $

#include "mulAlphabet.h"
#include "distribution.h"
#include "errorMsg.h"
#include <iostream>
#include "logFile.h"


mulAlphabet::mulAlphabet(const alphabet* baseAlphabet, int mulFactor) :
_baseAlphabet(baseAlphabet->clone()),
_mulFactor(mulFactor),
_size(baseAlphabet->size() * mulFactor)
{}

mulAlphabet::mulAlphabet(const mulAlphabet& other) :
_baseAlphabet(other._baseAlphabet->clone()),
_mulFactor(other._mulFactor),
_size(other._size)
{}

mulAlphabet::~mulAlphabet()
{
	if (_baseAlphabet) delete (_baseAlphabet);
}

mulAlphabet& mulAlphabet::operator=(const mulAlphabet &other)
{
	if (_baseAlphabet) delete (_baseAlphabet);
	_baseAlphabet = other._baseAlphabet->clone();
	_mulFactor = other._mulFactor;
	_size = other._size;
	return (*this);
}

int mulAlphabet::unknown() const 
{
	return (convertFromBasedAlphaInt(_baseAlphabet->unknown()));
}

int mulAlphabet::gap() const
{
		return (convertFromBasedAlphaInt(_baseAlphabet->gap()));
}

int mulAlphabet::stringSize() const 
{
	return _baseAlphabet->stringSize();
}

bool mulAlphabet::isSpecific(const int id) const 
{
	if (id >= _size)
		return false;
	else
		return _baseAlphabet->isSpecific(convertToBasedAlphaInt(id));
}

/* The first _size characters should be first. The rest of the characters aren't multiplied.
For example, when using nucleotides as the based alphabet and _mulFactor = 2 :
0	A0
1	C0
2	G0
3	T0
4	A1
5	C1
6	G1
7	T1
8	A
9	C
10	G
11	T
12	U
13	R
14	Y
15	K
16	M
17	S
18	W
19	B
20	D
21	H
22	V
23	N
-1	-
*/

string mulAlphabet::fromInt(const int id) const 
{
	// category and categoryName are for debug purpose
	int category(_mulFactor);
	if (id>=0)
		category = min(id / _baseAlphabet->size() , _mulFactor) ; 
	string categoryName("");
	categoryName = int2string(category);
	int inCategoryId = convertToBasedAlphaInt(id);
	return (_baseAlphabet->fromInt(inCategoryId) + categoryName);
}

int mulAlphabet::convertFromBasedAlphaInt(int id) const
{
	if (id < 0)
		return (id);

	return (id + _size);
}

int mulAlphabet::fromChar(const string& str, const int pos) const 
{
	int id = _baseAlphabet->fromChar(str,pos);
	return (convertFromBasedAlphaInt(id));
}


vector<int> mulAlphabet::fromString(const string &str) const 
{
	vector<int> result = _baseAlphabet->fromString(str);
	vector<int>::iterator itr = result.begin();
	for (; itr != result.end(); ++itr)
		*itr = convertFromBasedAlphaInt(*itr);
	
	return (result);
}


int mulAlphabet::convertToBasedAlphaInt(int id) const
{
	if (id<0)
		return (id);
	if (id >= _size)
		return (id - _size);

	return (id % _baseAlphabet->size());
}



int mulAlphabet::relations(const int charInSeq, const int charToCheck) const
{ 
	int baseAlphabetSize = _baseAlphabet->size();
	int categoryInSeq(_mulFactor);
	if (charInSeq>=0)
		categoryInSeq = min(charInSeq/baseAlphabetSize , _mulFactor);
	
	int categoryToCheck(_mulFactor);
	if (charToCheck>=0)
		categoryToCheck = min(charToCheck/baseAlphabetSize , _mulFactor);
	
	if (categoryToCheck == _mulFactor)
		LOG(4,<<"mulAlphabet::relations charToCheck should belong to category < _mulFactor = " << _mulFactor << endl);

	if ((categoryInSeq == categoryToCheck) || (categoryInSeq == _mulFactor))
		return _baseAlphabet->relations(convertToBasedAlphaInt(charInSeq),convertToBasedAlphaInt(charToCheck));
	
	return 0;
}


int mulAlphabet::compareCategories(int charA, int charB) const
{
	int baseAlphabetSize = _baseAlphabet->size();
	int categoryA(_mulFactor);
	if (categoryA>=0)
		categoryA = min(charA/baseAlphabetSize,_mulFactor);

	int categoryB(_mulFactor);
	if (categoryB>=0)
		categoryB = min(charB/baseAlphabetSize,_mulFactor);

	if (categoryA<categoryB)
		return 1;
	else if (categoryB<categoryA)
		return -1;
	return (0);
}
