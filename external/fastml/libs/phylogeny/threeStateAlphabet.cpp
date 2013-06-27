#include "threeStateAlphabet.h"

threeStateAlphabet::threeStateAlphabet() {}

int threeStateAlphabet::fromChar(const char s) const{
	switch (s) {
		case '0': return 0; break;	
		case '1': return 1; break;	
		case '2': return 2; break;
		default:
			vector<string> err;
			err.push_back(" The threeStateAlphabet sequences contained the character: ");
			err[0]+=s;
			err.push_back(" threeStateAlphabet was not one of the following: ");
			err.push_back(" 0, 1, 2");
			errorMsg::reportError(err);
	}// end of switch
	return -99; // never suppose to be here.	
}// end of function

vector<int> threeStateAlphabet::fromString(const string &str) const {
	vector<int> vec;
	for (int i=0;i<str.size();i++)
		vec.push_back(fromChar(str[i]));
	return vec;
}

string threeStateAlphabet::fromInt(const int in_id) const{
	char res = 0;
	switch (in_id) {
		case 0 : res = '0'  ; break;
		case 1 : res = '1'  ; break;
		case 2 : res = '2'  ; break;
		default:
			vector<string> err;
			err.push_back("unable to print threeState_id. threeState_id was not one of the following: ");
			err.push_back("0,1,2");
			errorMsg::reportError(err);
	}//end of switch
	string vRes;
	vRes.append(1,res);
	return vRes;
}// end of function

// There are no relations here.
int threeStateAlphabet::relations(const int charInSeq, const int charToCheck) const{
	if (charInSeq == charToCheck)
		return 1;
	return 0;
}

int threeStateAlphabet::fromChar(const string& str, const int pos) const{
	return fromChar(str[pos]);
}




