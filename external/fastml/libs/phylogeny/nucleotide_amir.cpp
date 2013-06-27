// $Id: nucleotide.cpp 962 2006-11-07 15:13:34Z privmane $

#include "nucleotide.h"
#include "errorMsg.h"


/*nucleotide::nucleotide() {
	_relation.resize(4);
	for (int i=0; i < _relation.size(); ++i) {
		_relation[i].resize(16);
	}
	for (int s=0;s<4;++s) {
		for (int t=0;t<16;++t){
			_relation[s][t]=relationsInternal(s,t);
		}
	}
}
*/
nucleotide::nucleotide() {
	_relation.resize(5);
	for (int i=0; i < _relation.size(); ++i) {
		_relation[i].resize(17);
	}
	for (int s=0;s<5;++s) {
		for (int t=0;t<17;++t){
			_relation[s][t]=relationsInternal(s,t);
		}
	}
}
int nucleotide::fromChar(const string& str, const int pos) const {
	return fromChar(str[pos]);
}

vector<int> nucleotide::fromString(const string &str) const {
	vector<int> vec;
	for (int i=0;i<str.size();i++)
	  vec.push_back(fromChar(str[i]));
	return vec;
}

int nucleotide::fromChar(const char s) const {
	switch (s) {
		case 'A' : case'a' : return 0 ; break;//	A = adenine
		case 'C' : case'c' : return 1 ; break;//	C = cytosine
		case 'G' : case'g' : return 2 ; break;//	G = guanine
		case 'T' : case't' : return 3 ; break;//	T = thymine
		case '-' : case'_' : return 4 ; break; //	U = uracil
		//case 'U' : case'u' : return 4 ; break; //	U = uracil
		case 'R' : case'r' : return 5 ; break;//	R = purine	(same as [GA])
		case 'Y' : case'y' : return 6 ; break;//	Y = pyrimidine	(same as [TC]) 
		case 'K' : case'k' : return 7 ; break;//	K = keto	(same as [GT])
		case 'M' : case'm' : return 8 ; break;//	M = amino	(same as [AC])
		case 'S' : case's' : return 9 ; break;//	S = strong 	(same as [GC])
		case 'W' : case'w' : return 10; break;//	W = weak   	(same as [AT])
		case 'B' : case'b' : return 11; break;//	B =        	(same as [GTC])
		case 'D' : case'd' : return 12; break;//	D =        	(same as [GAT])
		case 'H' : case'h' : return 13; break;//	H =        	(same as [ACT])
		case 'V' : case'v' : return 14; break;//	V =        	(same as [GCA])
		case 'N' : case'n' : return 15; break;//	N = any 	(same as [ACGT])
		case '?' : case'*' : return 15; break;
		//case '-' : case'_' : return -1; break;
		case 'U' : case'u' : return 16; break;
		case 'x' : case'X' : return 15; break;
		case '.' : return -3; break;		// . is used in some sequence files as the character just in the line above...
		default:
		vector<string> err;
		err.push_back(" The nucleotide sequences contained the character: ");
		err[0]+=s;
		err.push_back(" The nucleotide was not one of the following: ");
		err.push_back("A, C, G, T, X, -, ?");
		err.push_back("a, c, g, t, x, _, *");
		errorMsg::reportError(err); 
	}
	return -99;
}

string nucleotide::fromInt(const int id) const {
	char x= fromIntInternal(id);
	string res;
	res.append(1,x);
	return res;
}

char nucleotide::fromIntInternal(const int in_id)  const {
	switch (in_id) {
		case 0 : return 'A'  ; break;
		case 1 : return 'C'  ; break;
		case 2 : return 'G'  ; break;
		case 3 : return 'T'  ; break;
		//case -1: return '-'  ; break;
		//case 16 : return '-'  ; break;
		//case 4  : return 'U'; break;
		case 16 : return 'U'  ; break;
		case 4  : return '-'; break;
		case 5  : return 'R'; break;
		case 6  : return 'Y'; break;
		case 7  : return 'K'; break;
		case 8  : return 'M'; break;
		case 9  : return 'S'; break;
		case 10 : return 'W'; break;
		case 11 : return 'B'; break;
		case 12 : return 'D'; break;
		case 13 : return 'H'; break;
		case 14 : return 'V'; break;
		case 15 : return 'N'; break;
		default:
			vector<string> err;
			err.push_back(" unable to print nucleotide. nucleotide was not one of the following: ");
			err.push_back("A, C, G, T, -, ?");
			err.push_back("a, c, g, t, _, *");
			errorMsg::reportError(err); // make the program quit
	}//end of switch
	return '!' ;		// for the lousy compiler
}

int nucleotide::relationsInternal(const int ctc,const int charInSeq
					   ) const{ //ctc=charToCheck
  switch (charInSeq){
    case 0 : if (ctc==0) return 1 ; break;//	A = adenine
    case 1 : if (ctc==1) return 1 ; break;//	C = cytosine
    case 2 : if (ctc==2) return 1 ; break;//	G = guanine
    case 3 : if (ctc==3) return 1 ; break;//	T = thymine
    case 4 : if (ctc==4) return 1 ; break; //	U = uracil
    case 5 : if (ctc==2||ctc==0) return 1 ; break;//	R = purine	(same as [GA])
    case 6 : if (ctc==3||ctc==1) return 1 ; break;//	Y = pyrimidine	(same as [TC]) 
    case 7 : if (ctc==2||ctc==3) return 1 ; break;//	K = keto	(same as [GT])
    case 8 : if (ctc==0||ctc==1) return 1 ; break;//	M = amino	(same as [AC])
    case 9 : if (ctc==2||ctc==1) return 1 ; break;//	S =       	(same as [GC])
    case 10: if (ctc==0||ctc==3) return 1 ; break;//	W =        	(same as [AT])
    case 11: if (ctc==2||ctc==3||ctc==1) return 1 ; break;//	B =        	(same as [GTC])
    case 12: if (ctc==2||ctc==0||ctc==3) return 1 ; break;//	D =        	(same as [GAT])
    case 13: if (ctc==0||ctc==1||ctc==3) return 1 ; break;//	H =        	(same as [ACT])
    case 14: if (ctc==2||ctc==1||ctc==0) return 1 ; break;//	V =        	(same as [GCA])
    case 15: if (ctc==0||ctc==1||ctc==2||ctc==3) return 1 ; break;//	N = any 	(same as [ACGT])
	case 16 : if (ctc==16) return 1 ; break;//	gap  
  };
  return 0;
};

