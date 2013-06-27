// 	$Id: indel.cpp 962 2006-11-07 15:13:34Z privmane $	

#include "indel.h"

indel::indel() {}

int indel::fromChar(const char s) const{
	switch (s) {
		case 'x' : case'X' : return 0; break;	
		case '-' : case'_' : return 1; break;
		default:
			vector<string> err;
			err.push_back(" The indel sequences contained the character: ");
			err[0]+=s;
			err.push_back(" Indel was not one of the following: ");
			err.push_back(" -, X");
			err.push_back(" _, x");
			errorMsg::reportError(err);
	}// end of switch
	return -99; // never suppose to be here.	
}// end of function

vector<int> indel::fromString(const string &str) const {
	vector<int> vec;
	for (int i=0;i<str.size();i++)
	  vec.push_back(fromChar(str[i]));
	return vec;
}

string indel::fromInt(const int in_id) const{
  char res = 0;
	switch (in_id) {
		case 0 : res = 'X'  ; break;
		case 1 : res = '-'  ; break;
		default:
		vector<string> err;
		err.push_back("unable to print indel_id. indel_id was not one of the following: ");
		err.push_back("X, -");
		err.push_back("x, _");
		errorMsg::reportError(err);
		}//end of switch
	string vRes;
	vRes.append(1,res);
	return vRes;
}// end of function

// There are no relations here.
int indel::relations(const int charInSeq, const int charToCheck) const{
	if (charInSeq == charToCheck)
		return 1;
	return 0;
}

int indel::fromChar(const string& str, const int pos) const{
	return fromChar(str[pos]);
}


