// $Id: amino.cpp 2414 2007-10-08 14:34:42Z adist $

#include "amino.h"

//VVint amino::_relation;

amino::amino() {
	_relation.resize(24); 	// relation should realy be an allocted, two dimentional array, not a vector. 
	for (int i=0; i < _relation.size(); ++i) { // this implementation would be much faster.  with some c-tricks, this checkup could be done with one access only.
		_relation[i].resize(20);
	}

	for (int k=-2;k<=21;++k){
		for (int j=0;j<20;++j){
			_relation[k+2][j]=relations_internal(k,j);
		}
	}
}

int amino::fromChar(const char s) const{
	switch (s) {
	case 'A' : case'a' : return 0 ; break;
	case 'R' : case'r' : return 1 ; break;
	case 'N' : case'n' : return 2 ; break;
	case 'D' : case'd' : return 3 ; break;
	case 'C' : case'c' : return 4 ; break;
	case 'Q' : case'q' : return 5 ; break;
	case 'E' : case'e' : return 6 ; break;
	case 'G' : case'g' : return 7 ; break;
	case 'H' : case'h' : return 8 ; break;
	case 'I' : case'i' : return 9 ; break;
	case 'L' : case'l' : return 10; break;
	case 'K' : case'k' : return 11; break;
	case 'M' : case'm' : return 12; break;
	case 'F' : case'f' : return 13; break;
	case 'P' : case'p' : return 14; break;
	case 'S' : case's' : return 15; break;
	case 'T' : case't' : return 16; break;
	case 'W' : case'w' : return 17; break;
	case 'Y' : case'y' : return 18; break;
	case 'V' : case'v' : return 19; break;
	case 'B' : case'b' : return 20 ; break; // aspartate(D) or asparagine(N)
	case 'Z' : case'z' : return 21 ; break; // glutamate (E) or glutamine(Q)
	case '-' : case'_' : return -1; break;
	case '?' : case'*' : return -2; break;
	case 'x' : case'X' : return -2; break;
	case '.' : return -3; break;
	default:
	  vector<string> err;
	  err.push_back(" The amino-acid sequences contained the character: ");
	  err[0]+=s;
	  err.push_back(" Amino acid was not one of the following: ");
	  err.push_back(" A, B, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, X, Z, -, ?");
	  err.push_back(" a, b, r, n, d, c, q, e, g, h, i, l, k, m, f, p, s, t, w, y, v, x, z, _, *");
	  errorMsg::reportError(err);
	}// end of switch
	return -99; // never suppose to be here.	
}// end of function

vector<int> amino::fromString(const string &str) const {
	vector<int> vec;
	for (int i=0;i<str.size();i++)
	  vec.push_back(fromChar(str[i]));
	return vec;
}

string amino::fromInt(const int in_id) const{
  char res = 0;
	switch (in_id) {
		case 0 : res = 'A'  ; break;
		case 1 : res = 'R'  ; break;
		case 2 : res = 'N'  ; break;
		case 3 : res = 'D'  ; break;
		case 4 : res = 'C'  ; break;
		case 5 : res = 'Q'  ; break;
		case 6 : res = 'E'  ; break;
		case 7 : res = 'G'  ; break;
		case 8 : res = 'H'  ; break;
		case 9 : res = 'I'  ; break;
		case 10: res = 'L'  ; break;
		case 11: res = 'K'  ; break;
		case 12: res = 'M'  ; break;
		case 13: res = 'F'  ; break;
		case 14: res = 'P'  ; break;
		case 15: res = 'S'  ; break;
		case 16: res = 'T'  ; break;
		case 17: res = 'W'  ; break;
		case 18: res = 'Y'  ; break;
		case 19: res = 'V'  ; break;
		case 20: res = 'B'  ; break;
		case 21: res = 'Z'  ; break;
		case -1: res = '-'  ; break;
		case -2: res = 'X'  ; break;
		case -3: res = '.'  ; break;
		default:
		vector<string> err;
		err.push_back(" unable to print amino ac_id. amino ac_id was not one of the following: ");
		err.push_back("A, B, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, Z, -, ?");
		err.push_back("a, b, r, n, d, c, q, e, g, h, i, l, k, m, f, p, s, t, w, y, v, z, _, *");
		errorMsg::reportError(err);
		}//end of switch
	string vRes;
	vRes.append(1,res);
	return vRes;
}// end of function

int amino::relations(const int charInSeq, const int charToCheck) const{
	if (charInSeq == -1) {
		errorMsg::reportError("gaps in the sequences. Either change gaps to ? or remove gap positions");
	}
	return _relation[charInSeq+2][charToCheck];// <-MATAN, HERE YOU SWITHCED THE ORDER...
}

int amino::fromChar(const string& str, const int pos) const{
	return fromChar(str[pos]);
}

int amino::relations_internal(const int charInSeq, const int charToCheck) const{
	if (charInSeq == charToCheck) return 1;
	else if (charInSeq == fromChar('?')) return 1;
	else if ((charInSeq == fromChar('B')) && 
		 ((charToCheck == fromChar('N')) || 
		  (charToCheck == fromChar('D')))) return 1; // B is either N or D
	else if ((charInSeq == fromChar('Z')) && 
		 ((charToCheck == fromChar('Q')) ||
		  (charToCheck == fromChar('E')))) return 1; // Z is either E or Q
	return 0;
}


vector<int> aminoUtility::codonOf(const int a, codon &cod){
	vector<int> codons;
	amino amin;
	string strAmino=amin.fromInt(a);
	map <string, string> genCode=cod.geneticCode();
	map <string, string>::iterator it=genCode.begin();
	int tmp2=genCode.size();
	while (it!=genCode.end()){
		string tmp=(*it).second;
		if ((*it).second==strAmino){
			string strCodon=(*it).first;
			int c=cod.fromChar(strCodon,0);
			codons.push_back(c);		
		}
		it++;
	}
	if (codons.empty()){
		cout<<tmp2<<" amino is  = "<<a<<endl;
		errorMsg::reportError("error in function aminoUtility::codonOf: no codon found for amino acid");
	}
	return codons;
}
