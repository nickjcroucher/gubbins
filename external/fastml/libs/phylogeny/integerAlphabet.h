#ifndef ___INTEGER_ALPH
#define ___INTEGER_ALPH

#include "alphabet.h"
#include "errorMsg.h"


class integerAlphabet : public alphabet {
public:
	explicit integerAlphabet(int size): _size(size){}; 
	virtual ~integerAlphabet() {}
	virtual alphabet* clone() const { return new integerAlphabet(*this); }
	int unknown() const  {return -2;}
	int gap() const  {errorMsg::reportError("The method integerAlphabet::gap() is used"); return -1;} 
	int size() const {return _size;}
	int stringSize() const; // one letter code.
	int relations(const int charInSeq, const int charToCheck) const;
	int fromChar(const string& str, const int pos) const;
	int fromChar(const char s) const;
	string fromInt(const int in_id) const;
	vector<int> fromString(const string& str) const;
	bool isSpecific(const int id) const {return true;}

private:
	int _size;

};

#endif
