// $Id: someUtil.cpp 6055 2009-04-03 21:19:38Z rubi $

#include "someUtil.h"
#include "errorMsg.h"
#include "talRandom.h"
#include <cmath>
#include <ctime>
#include <iterator>
#include <algorithm>
#include <string>
#include <cctype>
#include <cassert>
using namespace std;

// for the _mkdir call
#if defined(WIN32) || defined(SunOS) || defined(solaris)
  #include <direct.h>
#else
  #include <sys/file.h>
  #include <dirent.h>
//  #include <io.h>
#endif

//swap between the 4 variables such that the first becomes the second, second becomes the third and third becomes the fourth.
//used in functoin mnbrack below.
void shift3(MDOUBLE &a, MDOUBLE &b, MDOUBLE &c, const MDOUBLE d) {
	a=b;
	b=c;
	c=d;
}

MDOUBLE computeAverage(const vector<int>& vec) {
	MDOUBLE sum=0.0;
	for (int i=0; i < vec.size(); ++i) {
		sum+=static_cast<MDOUBLE>(vec[i]);
	}
	return sum/static_cast<MDOUBLE>(vec.size());
}

// X ~ Poisson(lamda) -->	P(X=k) = ((lamda^k)/k!) * e^(-lamda)
// It isn't smart to first calculate factorial(k) because the size of long int limits this calculation to k<=13
MDOUBLE copmutePoissonProbability(const int& k, const long double& lamda)
{
	assert(k>=0);
	long double tmp = pow(lamda,k); // tmp = (lamda^k)/k!
	
	for (int i=2; i<=k; ++i)
		tmp/=i;

	return (tmp * exp(-lamda));
}


MDOUBLE computeAverage(const vector<MDOUBLE>& vec) {
	MDOUBLE sum=0.0;
	for (int i=0; i < vec.size(); ++i) sum+=vec[i];
	return sum/static_cast<MDOUBLE>(vec.size());
}

MDOUBLE computeStd(const vector<int>& vec) {// page 60, Sokal and Rohlf
	MDOUBLE sum=0.0;
	MDOUBLE sumSqr=0.0;
	MDOUBLE vecSize = static_cast<MDOUBLE>(vec.size());
	for (int i=0; i < vec.size(); ++i) {
		sum+=static_cast<MDOUBLE>(vec[i]);
		sumSqr+=(static_cast<MDOUBLE>(vec[i])*static_cast<MDOUBLE>(vec[i]));
	}
	MDOUBLE res= sumSqr-(sum*sum/vecSize);
	res /= (vecSize-1.0);
	res = sqrt(res);
	return res;
}

MDOUBLE computeStd(const vector<MDOUBLE>& vec) {// page 60, Sokal and Rohlf
	MDOUBLE sum=0.0;
	MDOUBLE sumSqr=0.0;
	MDOUBLE vecSize = static_cast<MDOUBLE>(vec.size());
	for (int i=0; i < vec.size(); ++i) {
		sum+=vec[i];
		sumSqr+=(vec[i]*vec[i]);
	}
	MDOUBLE res= sumSqr-(sum*sum/vecSize);
	res /= (vecSize-1.0);
	res = sqrt(res);
	return res;
}

void computeRelativeFreqsFollowingOneChanged(MDOUBLE  newValFreq, int indexNewFreq,Vdouble &freqs){
	MDOUBLE proportionAfterOptimization = 1.0 - newValFreq;
	MDOUBLE proportionBeforeOptimization = 1.0 - freqs[indexNewFreq];
	MDOUBLE sum = 0.0;
	for (int i=0; i<freqs.size(); ++i) {
		if (i==indexNewFreq){
			freqs[i] = newValFreq;
		}
		else {
			freqs[i] = proportionAfterOptimization*freqs[i]/proportionBeforeOptimization;
		}
		sum+=freqs[i];
	}
	if (!DEQUAL(sum,1.0)) {
		errorMsg::reportError("Error in computeRelativeFreqsFollowingOneChanged, sum not equal to 1");
	}
}


char mytolower(char in){return tolower(in);}
char mytoupper(char in){return toupper(in);}

void toLower(string& str) {
	transform (str.begin(), str.end(), str.begin(), mytolower);
}
void toUpper(string& str) {
	transform (str.begin(), str.end(), str.begin(), mytoupper);
}
bool allowCharSet(const string& allowableChars, const string& string2check) {
// this function check if all the character in string2check are made of characters from allowableChars
	for (int i=0; i < string2check.size(); ++i) {
		// now checking for string2check[i]
		int j;
		for (j=0; j < allowableChars.size(); ++j) {
			if (string2check[i] == allowableChars[j]) {
				break;
			}
		}
		if (j==allowableChars.size()) return false;
	}
	return true;
}

bool isCharInString(const string& stringToCheck, const char charToCheck) {
	for (int i=0; i < stringToCheck.size(); ++i ) {
		if (stringToCheck[i] == charToCheck) return true;
	}
	return false;
}

string double2string(const double x, const int lenght){
	
	// first getting the integer part:
	//Itay: fixing bug regarding negative floats 
	double x_abs = fabs(x);
	int theIntegerPart = static_cast<int>(x_abs);
	double theRemainingPart = fabs(x_abs-theIntegerPart);
	int integerRepresentingTheRemainingPart = static_cast<int>(theRemainingPart*pow(10.0,lenght));
	string part1 = int2string(theIntegerPart);
	string part2 = int2string(integerRepresentingTheRemainingPart);
	while (part2.length()<lenght){
		part2.insert(0, "0");
	}

	string result("");
	if (x < 0.0)
		result += "-";
	result += part1;
	result += ".";
	result += part2;

	// removing 0 from the end
	int i = result.length()-1;
	while (result[i]!='.' && i>0 && result[i]=='0'){
		result.erase(i);
		i--;
	}
	
	// removing "." if this is the last character in the string.
	if (result[result.length()-1]=='.')
	result.erase(result.length()-1);

	return result;
}

string int2string(const int num) {
// the input to this program is say 56
// the output is the string "56"
// this version of int2string is more portable 
// than sprintf like functions from c;
// or sstream of stl.
	if (num == 0) return "0";
	string res;
	int i = abs(num);
	
	
	int leftover;
	char k;
	while (i) {
		leftover = i%10;
		k = '0'+leftover;
		res = k+res;
		i/=10;
	}
	if (num<0) res = "-" + res;
	return res;
};

void printTime(ostream& out) {
	time_t ltime;
	time( &ltime );
	out<<"# the date is "<< ctime( &ltime )<<endl;
}

MDOUBLE string2double(const string& inString) {

	if (allowCharSet("0123456789.eE+-",inString) == false) {
		errorMsg::reportError(" error in function string2double ");
	}
	
	// first decide if the format is like 0.00343 (regularFormat) or
	// if it is in the form of 0.34e-006 for example

	bool regularFormat = true;
	int i;
	for (i=0; i < inString.size(); ++i) {
		if ((inString[i] == 'e' ) || (inString[i] == 'E' )) {
			regularFormat = false; 
			break;
		}
	}

	if (regularFormat) {
			MDOUBLE dDistance = atof(inString.c_str());
			return dDistance;
	}
	else {
		string b4TheExp;
		bool plusAfterTheExp = true;
		string afterTheExp;

		// b4 the exp
		for (i=0; i < inString.size(); ++i) {
			if (inString[i] != 'e' ) {
				b4TheExp += inString[i];
			}
			else break;
		}
		++i; //now standing after the exp;
		if (inString[i] == '-' ) {
			plusAfterTheExp = false;
			++i;
		}
		else if (inString[i] == '+' ) {
			plusAfterTheExp = true;
			++i;
		}
		else plusAfterTheExp = true; // the number is like 0.34e43

		for (; i < inString.size(); ++i) {
			afterTheExp += inString[i];
		}

		MDOUBLE res = 0.0;
		MDOUBLE dDistance = atof(b4TheExp.c_str());
		int exponentialFactor = atoi(afterTheExp.c_str());
		if (plusAfterTheExp) res = dDistance * pow(10.0,exponentialFactor);
		else res = dDistance * pow(10.0,-exponentialFactor);

		return res;
	}

	
}


bool checkThatFileExist(const string& fileName) {
	ifstream file1(fileName.c_str());
	if (file1==NULL) return false;
	file1.close();
	return true;
}

void putFileIntoVectorStringArray(istream &infile,vector<string> &inseqFile){
	inseqFile.clear();
	string tmp1;
	while (getline(infile,tmp1, '\n' ) ) {
		if (tmp1.size() > 15000) {
			vector<string> err;
			err.push_back("Unable to read file. It is required that each line is no longer than");
			err.push_back("15000 characters. ");
			errorMsg::reportError(err,1); 
		}
		if (tmp1[tmp1.size()-1]=='\r') {// in case we are reading a dos file 
			tmp1.erase(tmp1.size()-1);
		}// remove the traling carrige-return
		inseqFile.push_back(tmp1);
	}
}

bool fromStringIterToInt(string::const_iterator & it, // ref must be here
						const string::const_iterator endOfString,
						int& res) {// the ref is so that we can use the it after the func.
	while (it != endOfString) {
		if ((*it == ' ') || (*it == '\t')) ++it;else break; // skeeping white spaces.
	}
	if (it != endOfString) {
		if (isdigit(*it) || (*it == '-')){
			int k = atoi(&*it);
			if (*it == '-') ++it;
			for (int numDig = abs(k); numDig>0; numDig/=10) ++it;
			res = k;
			return true;
		}
		else return false; //unable to read int From String
	}
	return false; //unable to read int From String
	
}

string* searchStringInFile(const string& string2find,
						   const int index,
						   const string& inFileName) {
	ifstream f;
	f.open(inFileName.c_str());
	if (f==NULL) {
		string tmp = "Unable to open file name: "+inFileName+" in function searchStringInFile"; 
		errorMsg::reportError(tmp);
	}

	string numm = int2string(index);
	string realString2find = string2find+numm;

	istream_iterator<string> is_string(f);
	istream_iterator<string> end_of_stream;

	is_string = find(is_string,end_of_stream,realString2find);
	if(is_string == end_of_stream) {f.close();return NULL;}
	else {
		is_string++;
		if(is_string == end_of_stream) {f.close();return NULL;};
		string* s = new string(*is_string);
		f.close();
		return s;
	}
	f.close();
	return NULL;
}
string* searchStringInFile(const string& string2find,
						   const string& inFileName) {// return the string that is AFTER the string to search.
	ifstream f;
	f.open(inFileName.c_str());
	if (f==NULL) {
		string tmp = "Unable to open file name: "+inFileName+" in function searchStringInFile"; 
		errorMsg::reportError(tmp);
	}
	string realString2find = string2find;

	istream_iterator<string> is_string(f);
	istream_iterator<string> end_of_stream;

	is_string = find(is_string,end_of_stream,realString2find);
	if(is_string == end_of_stream) {f.close();return NULL;}
	else {
		is_string++;
		if(is_string == end_of_stream) {f.close();return NULL;};
		string* s = new string(*is_string);
		f.close();
		return s;
	}
	f.close();
	return NULL;
}
bool doesWordExistInFile(const string& string2find,const string& inFileName) {
	ifstream f;
	f.open(inFileName.c_str());
	if (f==NULL) {
		string tmp = "Unable to open file name: "+inFileName+" in function searchStringInFile"; 
		errorMsg::reportError(tmp);
	}

	istream_iterator<string> is_string(f);
	istream_iterator<string> end_of_stream;

	is_string = find(is_string,end_of_stream,string2find);
	if(is_string == end_of_stream) return false;
	else return true;
}

string takeCharOutOfString(const string& charsToTakeOut, const string& fromString) {
	string finalString;
	for (int i=0; i<fromString.size(); ++i) {
		bool goodChar = true;
		for (int j=0; j < charsToTakeOut.size(); ++j) {
			if (fromString[i]== charsToTakeOut[j]) goodChar = false;
		}
		if (goodChar) finalString+=fromString[i];
	}
	return finalString;
}

bool DEQUAL(const MDOUBLE x1, const MDOUBLE x2, MDOUBLE epsilon/*1.192092896e-07F*/) {
	return (fabs(x1-x2)<epsilon);
}

bool DBIG_EQUAL(const MDOUBLE x1, const MDOUBLE x2, MDOUBLE epsilon/*1.192092896e-07F*/){
	return ((x1 > x2) || DEQUAL(x1, x2,epsilon));
}


bool DSMALL_EQUAL(const MDOUBLE x1, const MDOUBLE x2, MDOUBLE epsilon/*1.192092896e-07F*/){ 
	return ((x1 < x2) || DEQUAL(x1, x2,epsilon));
}

void createDir(const string & curDir, const string & dirName){// COPYRIGHT OF ITAY MAYROSE.
	string newDir;
	if (curDir == "")
		newDir = dirName;
	else
		newDir = curDir + string("/") + dirName;
#ifdef WIN32
	if( _mkdir(newDir.c_str()) == 0 ){
		LOG(5, << "Directory " <<newDir<<" was successfully created"<<endl);
    }else{
		if (errno == EEXIST) {
			LOG(5,<<"Directory already exist");
			return;
		} else {
		string err = "Problem creating directory " + newDir + " \n";
		LOG(5, << err << endl);
		errorMsg::reportError(err);
		}
	}
#else
	DIR * directory = opendir(newDir.c_str());
	if (directory == NULL) {
		string sysCall = "mkdir " + newDir;
        	system(sysCall.c_str());
	}
	else{
		string err = "Directory " + newDir + " already exists  \n";
		LOG(5, << err << endl);
		//errorMsg::reportError(err);
		
	}
#endif
}

//scale vecToScale so that its new average is AvgIn. return the scaling factor. 
MDOUBLE scaleVec(Vdouble& vecToScale, const MDOUBLE avgIn)
{
	int vecSize = vecToScale.size();
	MDOUBLE sum = 0;
	for (int x = 0; x<vecSize; ++x)
	{
		sum += vecToScale[x];
	}
	MDOUBLE avg = sum/vecSize;
	MDOUBLE scaleFactor = avgIn / avg;
	
	for (int i = 0; i<vecSize; ++i)
	{
		vecToScale[i] *= scaleFactor;
	}

	MDOUBLE newAvg = computeAverage(vecToScale);
	if (fabs(newAvg - avgIn) > 0.001)
		errorMsg::reportError(" problem - scalled average is not avgIn after scalling!!!");
	return scaleFactor;
}

//calculates the mean square error distance between 2 vectors:
MDOUBLE calcMSEDistBetweenVectors(const Vdouble& oneRatesVec, const Vdouble& otherRatesVec)
{
	MDOUBLE res = 0.0;
	if (oneRatesVec.size() != otherRatesVec.size())
		errorMsg::reportError("the two vectors to be compared are not the same size in function SimulateRates::calcDistBetweenRatesVectors()");
	
	for (int i=0; i<oneRatesVec.size(); ++i)
	{
		MDOUBLE diff = oneRatesVec[i] - otherRatesVec[i];
		res += diff * diff;
	}

	res /= oneRatesVec.size();
	return res;
}

//calculates the mean absolute deviations distance between 2 vectors:
MDOUBLE calcMADDistBetweenVectors(const Vdouble& oneRatesVec, const Vdouble& otherRatesVec)
{
	MDOUBLE res = 0.0;
	if (oneRatesVec.size() != otherRatesVec.size())
		errorMsg::reportError("the two vectors to be compared are not the same size in function SimulateRates::calcDistBetweenRatesVectors()");
	
	for (int i=0; i<oneRatesVec.size(); ++i)
	{
		MDOUBLE diff = oneRatesVec[i] - otherRatesVec[i];
		res += fabs(diff);
	}

	res /= oneRatesVec.size();
	return res;
}

MDOUBLE calcRelativeMADDistBetweenVectors(const Vdouble& trueValues, const Vdouble& inferredValues, const MDOUBLE threshhold/*0.0*/)
{
	MDOUBLE res = 0.0;
	if (inferredValues.size() != trueValues.size())
		errorMsg::reportError("the two vectors to be compared are not the same size in function SimulateRates::calcDistBetweenRatesVectors()");
	
	int counter = 0;
	for (int i=0; i<inferredValues.size(); ++i)
	{
		if (trueValues[i] < threshhold)
			continue;
		MDOUBLE diff = fabs(inferredValues[i] - trueValues[i]);
		res += (diff / trueValues[i]);
		++counter;
	}

	res /= counter;
	return res;
}

//calculates the relative mean square error distance between 2 vectors:
//The difference from a regualar MSE is that for each position the squared difference is devided by the true value
//if threshhold > 0: if trueValues[i] < threshhold then do not add the rse for this psition to the result
MDOUBLE calcRelativeMSEDistBetweenVectors(const Vdouble& trueValues, const Vdouble& inferredValues, const MDOUBLE threshhold/*0.0*/ )
{
	MDOUBLE res = 0.0;
	if (inferredValues.size() != trueValues.size())
		errorMsg::reportError("the two vectors to be compared are not the same size in function SimulateRates::calcDistBetweenRatesVectors()");
	
	int counter = 0;
	for (int i=0; i<inferredValues.size(); ++i)
	{
		if (trueValues[i] < threshhold)
			continue;
		MDOUBLE diff = inferredValues[i] - trueValues[i];
		res += diff * diff / trueValues[i];
		++counter;
	}

	res /= counter;
	return res;
}


MDOUBLE calcRankCorrelation(const Vdouble& oneRatesVec, const Vdouble& otherRatesVec)
{
	MDOUBLE res = 0.0;
	Vdouble orderVec1, orderVec2;
	MDOUBLE s_one = orderVec(oneRatesVec, orderVec1);
	MDOUBLE s_two = orderVec(otherRatesVec, orderVec2);
	int seqLength = oneRatesVec.size();
	MDOUBLE diff, sum_diff_sqr = 0;
	for (int i=0; i<seqLength; ++i)
	{
		diff = orderVec1[i] - orderVec2[i];
		sum_diff_sqr += pow(diff, 2);
	}
	MDOUBLE en3n = (seqLength * (pow(static_cast<double>(seqLength), 2.0) -1)); //n^3 -n
	MDOUBLE numerator = 1.0 - ((6/en3n) * (sum_diff_sqr + (s_one + s_two)/12.0));
	MDOUBLE denum = sqrt((1.0 - s_one/en3n) * (1.0 - s_two/en3n));
	res = numerator/ denum;
	return res;
}

ostream &operator<<(ostream &out, const Vdouble &v){
  for (int j=0;j<v.size();++j)
    out<< v[j]<<" ";
  out <<endl;
  return(out);
}

ostream &operator<<(ostream &out, const VVdouble &m){
  for (int i=0;i<m.size();++i)
    out<<m[i];
  out <<endl;
  return(out);
}

void mult(Vdouble& vec, const MDOUBLE factor){
  for(int i=0;i<vec.size();++i)
    vec[i]*=factor;
}

void mult(VVdouble& vec, const MDOUBLE factor){
  for(int i=0;i<vec.size();++i)
    mult(vec[i],factor);
}



////orderVec - determine the relative order of vecIn
////returns orderVecOut[i] is the rank of vecIn[i]
////note that in case of ties the rank will be the midrank of the tied group
//Vdouble orderVec(const Vdouble& vecIn)
//{
//	int vecSize = vecIn.size();
//	Vdouble orderVecOut(vecSize);
//	vector< vecElem<MDOUBLE> > sortVec(vecSize);
//	for (int x =0; x < vecSize ; ++x)
//	{
//		sortVec[x].setValue(vecIn[x]);
//		sortVec[x].setPlace(x);
//	}
//	sort(sortVec.begin(), sortVec.end());
//	
//	//check for ties and correct their rank
//	Vdouble rankVec(vecSize);
//	MDOUBLE rank;
//	for (int i=0; i < vecSize; )
//	{
//		if (sortVec[i].getValue() != sortVec[i+1].getValue())
//		{//no tie
//			rankVec[i] = i;
//			++i;
//		}
//		else
//		{//tie
//			int to =0;
//			for (to = i+1; (to<=vecSize) && (sortVec[i].getValue() == sortVec[to].getValue());++to)
//				;//check how far the tie goes
//			to--;
//			rank = 0.5*(to + i); 
//			for (int ji = i; ji<= to; ji++)
//			{
//				rankVec[ji] = rank;
//			}
//
//			i = to+1;
//		}
//	}
//	for (int j =0; j < vecSize; ++j) {
//		assert ((rankVec[j] >= 0) && (rankVec[j] < vecSize));
//		orderVecOut[sortVec[j].getPlace()] = rankVec[j]; 
//	}
//	return orderVecOut;
//}

//orderVec - determine the relative order of vecIn
//orderVecOut[i] is the rank of vecIn[i]
//note that in case of ties the rank will be the midrank of the tied group
//return sum of n^3 - n where n is the number of elements in each tied group - see spearman rank correlation
MDOUBLE orderVec(const vector<MDOUBLE>& vecIn, vector<MDOUBLE>& orderVecOut)
{
	int vecSize = vecIn.size();
	orderVecOut.resize(vecSize);
	vector< vecElem<MDOUBLE> > sortVec(vecSize);
	for (int x =0; x < vecSize ; ++x)
	{
		sortVec[x].setValue(vecIn[x]);
		sortVec[x].setPlace(x);
	}
	sort(sortVec.begin(), sortVec.end());
	//check for ties and correct their rank
	Vdouble rankVec(vecSize);
	MDOUBLE sumRankDif = 0; //sum(Fk^3 - Fk)

	MDOUBLE rank;
	for (int i=0; i < vecSize; )
	{
		if (sortVec[i].getValue() != sortVec[i+1].getValue())
		{//no tie
			rankVec[i] = i;
			++i;
		}
		else
		{//tie
			int to =0;
			for (to = i+1; (to<=vecSize) && (sortVec[i].getValue() == sortVec[to].getValue());++to)
				;//check how far the tie goes
			to--;
			rank = 0.5*(to + i); 
			for (int ji = i; ji<= to; ji++)
			{
				rankVec[ji] = rank;
			}

			int numTies = to - i +1; //number o fties in this group
			sumRankDif += numTies*numTies*numTies - numTies;
			i = to+1;
		}
	}
	
	for (int j =0; j < vecSize; ++j) {
	    assert ((rankVec[j] >= 0) && (rankVec[j] < vecSize));
	    orderVecOut[sortVec[j].getPlace()] = rankVec[j]; 
	}
	return sumRankDif;
}

void orderVec(const Vdouble& vecIn, vector< vecElem<MDOUBLE> >& orderVecOut)
{
	int vecSize = vecIn.size();
	orderVecOut.resize(vecSize);
	for (int x =0; x < vecSize ; ++x)
	{
		orderVecOut[x].setValue(vecIn[x]);
		orderVecOut[x].setPlace(x);
	}
	sort(orderVecOut.begin(), orderVecOut.end());
}


void splitString2(string str, string seperater, string &first, string &second)
{
     int i = (int)str.find(seperater); //find seperator
     if(i != -1)
     {
          int y = 0;
          if(!str.empty())
          {
               while(y != i)
               {
                    first += str[y++]; //creating first string
               }
               y = y+(int)seperater.length(); //jumping forward seperater length
               while(y != str.length())
               {
                    second += str[y++]; //creating second string
               }
               
          }
     }
     else
     {
          first = str;
          second = "NULL"; //if seperator is not there then second string == null 
     }
}


void splitString(const string& str,vector<string>& subStrs,const string& delimiter)
{
	// Skip delimiter at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiter,0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiter,lastPos);

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a subStr, add it to the vector.
		subStrs.push_back(str.substr(lastPos,pos - lastPos));
		// Skip delimiter.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiter,pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiter,lastPos);
	}
}

Vint getVintFromStr(const string& inStr)
{
	Vint res;
	vector<string> outStr;
	splitString(inStr, outStr, ",");
	for (int i = 0; i < outStr.size(); ++i)
	{
		int x = atoi(outStr[i].c_str());
		res.push_back(x);
	}
	return res;
}

string getStrFromVint(const Vint& inVec)
{
	string res("");
	for (int i = 0; i < inVec.size(); ++i)
	{
		if (i > 0)
			res += ",";
		res += int2string(inVec[i]);
	}
	return res;
}


/********************************************************************************************
*********************************************************************************************/
int fromIndex2gainIndex(const int i, const int gainCategories, const int lossCategories){
	int gainIndex;
	if(lossCategories<=gainCategories){			
		gainIndex = (int)floor((double)i/(lossCategories) );
	}
	else{
		gainIndex = i%(gainCategories);			
	}
	return gainIndex;
}

int fromIndex2lossIndex(const int i, const int gainCategories, const int lossCategories){
	int lossIndex;
	if(lossCategories<=gainCategories){			
		lossIndex = i%(lossCategories);
	}
	else{
		lossIndex = (int)floor((double)i/(gainCategories) );		
	}
	return lossIndex;
}

int giveRandomState(const int alphabetSize, const int beginningState, const VVdouble &changeProbabilities) 
{
	for (int loop = 0 ; loop < 100000 ; ++loop) 
	{
		MDOUBLE theRandNum = talRandom::giveRandomNumberBetweenZeroAndEntry(1.0);
		MDOUBLE sum = 0.0;
		for (int state = 0; state < alphabetSize; ++state) 
		{
			sum += changeProbabilities[beginningState][state];
			if (theRandNum < sum) {
				return state;
			}
		}
	}
	errorMsg::reportError("giveRandomState: could not give random character. The reason is unknown.");
	return 1;

}

int giveRandomState(const int alphabetSize, const Vdouble &frequencies)  {
	for (int loop =0 ;loop<100000 ;loop++) {
		MDOUBLE theRandNum = talRandom::giveRandomNumberBetweenZeroAndEntry(1.0);
		MDOUBLE sum = 0.0;
		for (int j=0; j < alphabetSize;++j) {
			sum+=frequencies[j];
			if (theRandNum<sum) return j;
		}
	} 
	errorMsg::reportError("giveRandomState: Could not give random character. The reason is probably that the frequencies do not sum to one.");
	return 1;
}

