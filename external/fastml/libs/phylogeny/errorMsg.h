// $Id: errorMsg.h 962 2006-11-07 15:13:34Z privmane $

// version 1.01
// last modified 1 Jan 2004

#ifndef ___ERROR_MSG_H
#define ___ERROR_MSG_H

#include <string>
#include <vector>
#include <iostream>

using namespace std;

// The error is always send to cerr. _errorOut is NULL, unless setErrorOstream is called.


class errorMsg {
public:
	static void reportError(const vector<string>& textToPrint, const int exitCode=1);
	static void reportError(const string& textToPrint, const int exitCode=1);
	static void setErrorOstream(ostream* errorOut) {_errorOut = errorOut;}
private:
	static ostream* _errorOut;
};

// example of how to output to a file called error.txt
// ofstream f("error.txt");
// errorMsg::setErrorOstream(&f);
// errorMsg::reportError("cheers");

#endif

