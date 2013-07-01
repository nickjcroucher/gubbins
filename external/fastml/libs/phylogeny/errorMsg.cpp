// $Id: errorMsg.cpp 6066 2009-04-14 19:11:10Z itaymay $

// version 1.01
// last modified 1 Jan 2004
#include "definitions.h"
#include <cassert>
#include "errorMsg.h"
#include "logFile.h"
#include <errno.h>
#include <string.h> //for strerror
#include <stdlib.h>  //for exit()

ostream *errorMsg::_errorOut= NULL;

void errorMsg::reportError(const vector<string>& textToPrint, const int exitCode) {
	for (int i =0 ; i < textToPrint.size() ; ++i) {
		LOG(1,<<textToPrint[i]<<endl);
		cerr<<textToPrint[i]<<endl;
		if (_errorOut != NULL && *_errorOut != cerr)  {
			(*_errorOut)<<textToPrint[i]<<endl;
		}
	}
	if (errno!=0){
	  LOG(1,<<"System Error: "<<strerror(errno)<<endl);
	  cerr<<"System Error: "<<strerror(errno)<<endl;
	}
	assert(0); // always stop here if in DEBUG mode.
	exit(exitCode);
}

void errorMsg::reportError(const string& textToPrint, const int exitCode) {
	LOG(1,<<endl<<textToPrint<<endl);
	cerr<<endl<<textToPrint<<endl;
	if (_errorOut != NULL && *_errorOut != cerr)  {
		(*_errorOut)<<textToPrint<<endl;
	}
	if (errno!=0){
	  LOG(1,<<"System Error: "<<strerror(errno)<<endl);
	  cerr<<"System Error: "<<strerror(errno)<<endl;
	}
	assert(0); // always stop here if in DEBUG mode.
	exit(exitCode);
}


