// $Id: someUtil.h 6055 2009-04-03 21:19:38Z rubi $

#ifndef ___SOME_UTIL_H
#define ___SOME_UTIL_H

#include "logFile.h"
#include "definitions.h"
#include "alphabet.h"
#include <string>
#include <iostream>
using namespace std;

//to be used for orderVec
template <class T>
class vecElem
{
public:
	vecElem();
	virtual ~vecElem() {};
	void setValue(const T val) {m_value = val;}
	T getValue() {return m_value;}
	void setPlace(const int place) {m_place = place;}
	int getPlace() {return m_place;}
	inline bool operator< (const vecElem& elemIn) const;
private:
	int m_place;
	T m_value;
};


template <class T>
vecElem< T >::vecElem()
{
	m_value = -1;
	m_place = -1;
}

//template <class T>
//vecElement< T >::~vecElement()
//{
//}
template <class T>
bool vecElem< T >::operator<(const vecElem& elemIn) const
{
	if (m_value == elemIn.m_value)
		return (m_place  < elemIn.m_place);
	else
		return (m_value < elemIn.m_value);
}



// STATISTICAL UTILITIES:

MDOUBLE computeAverage(const vector<int>& vec);
MDOUBLE computeAverage(const vector<MDOUBLE>& vec);
MDOUBLE computeStd(const vector<MDOUBLE>& vec);// page 60, Sokal and Rohlf
MDOUBLE computeStd(const vector<int>& vec);// page 60, Sokal and Rohlf
MDOUBLE copmutePoissonProbability(const int& k, const long double& lamda);
// re-computes a vector of frequencies after one value is changed: 
// all other values are set according to their relative value 
void computeRelativeFreqsFollowingOneChanged(MDOUBLE newValFreq, int indexNewFreq,Vdouble &freqs);//freqs is the old vector into which we write the new values

// SIMULATIONS:
int giveRandomState(const int alphabetSize, const int beginningState, const VVdouble &changeProbabilities);
int giveRandomState(const int alphabetSize, const Vdouble &frequencies);

// TIME UTILITIES
void printTime(ostream& out);

// TEXT UTILITIES
string int2string(const int i);
string double2string(const double x, int const howManyDigitsAfterTheDot=5);
MDOUBLE string2double(const string& inString);
bool allowCharSet(const string& allowableChars, const string& string2check);
bool isCharInString(const string& stringToCheck, const char charToCheck);
void putFileIntoVectorStringArray(istream &infile,vector<string> &inseqFile);

bool fromStringIterToInt(string::const_iterator & it,
						 const string::const_iterator endOfString,
						 int& res);

string takeCharOutOfString(const string& charsToTakeOut, const string& fromString);
void toLower(string& str);
void toUpper(string& str);
//splits the string to substr according to the given delimiter (parallel to split in perl)
void splitString(const string& str,vector<string>& subStrs,const string& delimiter);

//input: a list of INTs seperated by commas ("1,3,5") returns the int in the vector
Vint getVintFromStr(const string& str);
//return a list of INTs seperated by commas ("1,3,5") 
string getStrFromVint(const Vint& inVec);

// FILE UTILITIES
bool checkThatFileExist(const string& fileName); 
string* searchStringInFile(const string& string2find,
						   const int index,
						   const string& inFileName);
string* searchStringInFile(const string& string2find,
						   const string& inFileName);
bool doesWordExistInFile(const string& string2find,const string& inFileName);
void createDir(const string& curDir,const string& dirName);


//BIT UTILITIES
//void nextBit(bitset<64> &cur);

//ARITHMETIC UTILITIES
//DEQUAL: == UP TO EPSILON
//DBIG_EQUAL: >= UP TO EPSILON
//DSMALL_EQUAL: <= UP TO EPSILON
bool DEQUAL(const MDOUBLE x1, const MDOUBLE x2, const MDOUBLE epsilon = 1.192092896e-07F); // epsilon taken from WINDOW'S FILE FLOAT.H
bool DBIG_EQUAL(const MDOUBLE x1, const MDOUBLE x2, const MDOUBLE epsilon = 1.192092896e-07F); 
bool DSMALL_EQUAL(const MDOUBLE x1, const MDOUBLE x2, const MDOUBLE epsilon = 1.192092896e-07F); // {return ((x1 < x2) || DEQUAL(x1, x2));}

//swap between the 4 variables such that the first becomes the second, second becomes the third and third becomes the fourth.
//used in functoin mnbrack below.
void shift3(MDOUBLE &a, MDOUBLE &b, MDOUBLE &c, const MDOUBLE d);


// print vector and VVdoulbe util
ostream &operator<<(ostream &out, const Vdouble &v);
ostream &operator<<(ostream &out, const VVdouble &m);
void mult(Vdouble& vec, const MDOUBLE factor);
void mult(VVdouble& vec, const MDOUBLE factor);
//scale vecToScale so that its new average is AvgIn. return the scaling factor. 
MDOUBLE scaleVec(Vdouble& vecToScale, const MDOUBLE avgIn);
//determine the relative order of vecIn. The order vector is returned 
//ex: vecIn = [0.1 0.4 0.01 0.9 1.8] orderVecOut = [1 2 0 3 4] 
MDOUBLE orderVec(const vector<MDOUBLE>& vecIn, vector<MDOUBLE>& orderVecOut);
//in this version orderVecOut does not preserv the same order as vecIn. 
//orderVecOut[0] cotains the lowest score and it is stored in orderVecOut[0].getValue()
//The place in the original vector is stored in orderVecOut[0].getPlace()
void orderVec(const Vdouble& vecIn, vector< vecElem<MDOUBLE> >& orderVecOut);
//calculates the spearman rank correlation value
MDOUBLE calcRankCorrelation(const Vdouble& oneRatesVec, const Vdouble& otherRatesVec);
MDOUBLE calcRelativeMSEDistBetweenVectors(const Vdouble& trueValues, const Vdouble& inferredValues, const MDOUBLE threshhold = 0.0);
MDOUBLE calcMSEDistBetweenVectors(const Vdouble& trueValues, const Vdouble& inferredValues);
//MAD = mean absolute deviations distance 
MDOUBLE calcMADDistBetweenVectors(const Vdouble& oneRatesVec, const Vdouble& otherRatesVec);
MDOUBLE calcRelativeMADDistBetweenVectors(const Vdouble& trueValues, const Vdouble& inferredValues, const MDOUBLE threshhold = 0.0);


/* Will split a string into 2 by the given seperator
Example for usage:
      string a, b, c;
      a.assign("Hello world!");
      splitString2(a, " ", b, c);
      cout << "b = " << b << endl << "c = " << c << endl;
      //b == Hello
      //c == world!
*/
void splitString2(string str, string seperater, string &first, string &second);

int fromIndex2gainIndex(const int i, const int gainCategories, const int lossCategories);
int fromIndex2lossIndex(const int i, const int gainCategories, const int lossCategories);


 
#endif

