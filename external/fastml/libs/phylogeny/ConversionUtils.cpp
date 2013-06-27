#include "ConversionUtils.h"
#include "someUtil.h"
#include "errorMsg.h"

#include <cmath>

using namespace std;

void appendIntToString (string& ioString, const int inValue) {
	std::ostringstream o;
	o << ioString<< inValue;
	ioString = o.str();
}

string appendInt2string(const int x)
{
	string res;
	appendIntToString(res, x);
	return res;
}

string appendDouble2string(const double x, const int lenght){
	
	// first getting the integer part:
	int theIntegerPart = static_cast<int>(x);
	double theRemainingPart = fabs(x-theIntegerPart);
	int integerRepresentingTheRemainingPart = static_cast<int>(theRemainingPart*pow(10.0,lenght));
	string part1, part2; 
	appendIntToString(part1, theIntegerPart);
	appendIntToString(part2, integerRepresentingTheRemainingPart);
	while (part2.length()<lenght){
		part2.insert(0, "0");
	}

	string result = part1;
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

