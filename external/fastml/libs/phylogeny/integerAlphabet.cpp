#include "integerAlphabet.h"
#include "logFile.h"
#include "someUtil.h"
#include <cctype>

//return -99 if not succeeds.
int integerAlphabet::fromChar(const string& s, const int pos) const {
	if (s.size() <= (pos + stringSize()-1)) {
		string textToPrint("integerAlphabet::fromChar: Trying to read a character past the end of the string. ");
		LOG(1,<<textToPrint<<endl);
		return -99;
	}

	string s_sub=s.substr(pos,stringSize());
	int leftMostDigit(0);
	// find the left most digit. (s_sub can contain for example "0032" and so the left most digit is '3' and the number that should be returned is 32.
	for (leftMostDigit=0; leftMostDigit < s_sub.size(); ++leftMostDigit) {
		if (s_sub[leftMostDigit]!='0')
			break;
	}
	s_sub =s_sub.substr(leftMostDigit);

	return (atoi(s_sub.c_str()));
}

vector<int> integerAlphabet::fromString(const string &str) const {
	vector<int> vec;
	if (str.size()%stringSize()!=0) {
		errorMsg::reportError("error in integerAlphabet::fromString. String length should be a multiplication of stringSize");
	}
	for (int i=0;i<str.size();i+=stringSize())
	  vec.push_back(fromChar(str,i));
	return vec;
}


int integerAlphabet::stringSize() const {
	int countDigits(1);
	int wholeNum = _size/10;
	while (wholeNum > 0) {
		countDigits++;
		wholeNum /=10;
	}
	return (countDigits);
}


string integerAlphabet::fromInt(const int in_id) const{
	
	string res = int2string(in_id);
	while (res.size() <= stringSize()) {
	}
	return res;
}

// There are no relations here.
int integerAlphabet::relations(const int charInSeq, const int charToCheck) const{
	if (charInSeq == charToCheck)
		return 1;
	return 0;
}
