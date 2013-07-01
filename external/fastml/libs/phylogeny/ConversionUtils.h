//utility class that converts between data types
#ifndef ___ConversionUtils_h
#define ___ConversionUtils_h

#include <sstream>
#include <string>
#include "definitions.h"

using namespace std;

//a function that turns an integer to string 

void appendIntToString (string& ioString, const int inValue);
string appendDouble2string(const double x, int const howManyDigitsAfterTheDot=5);
string appendInt2string(const int x);


// Trims spaces at the left side of a string
static inline string trim_left(const string& str )
{
	int i=str.find_first_not_of(" \t");
	if(str.size()==0 || i >= str.size())
		return str;
	return str.substr( i ) ;
}
 

////
// Trims spaces at the right side of a string
static inline string trim_right(const string& str )
{
	int i=str.find_last_not_of(" \t");
	if(str.size()==0 || i >= str.size())
		return str;
	return str.substr(0, i + 1);
}
  
//// 
// Trims spaces at both sides of a string
static inline string trim(const string& str )
{
	return trim_left(trim_right(str));
}


#endif





