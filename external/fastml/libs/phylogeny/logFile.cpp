// $Id: logFile.cpp 962 2006-11-07 15:13:34Z privmane $

#include "logFile.h"
#include "errorMsg.h"

int myLog::_loglvl = 3;
ostream *myLog::_out= NULL;
bool myLog::_firstTime = true;

void myLog::setLog(const string logfilename, const int loglvl) {
	if (_out != NULL) myLog::endLog();
	if ((logfilename == "-")|| (logfilename == "")) {
		myLog::setLogOstream(&cout);
	} else {
		ofstream* outLF = new ofstream;
		if (_firstTime) {
		  outLF->open(logfilename.c_str());
		  _firstTime = false;
		}
		else
		  outLF->open(logfilename.c_str(), ofstream::out | ofstream::app); // append
		if (!outLF->is_open()) {
			errorMsg::reportError(string("Can't open for writing the log file ")+logfilename);
		}
		myLog::setLogOstream(outLF);
	}
	myLog::setLogLvl(loglvl);
	LOG(3,<<"START OF LOG FILE"<<endl);
}

void myLog::endLog(void){
	LOG(3,<<"END OF LOG FILE"<<endl);
        if (_out!=&cout && _out != NULL) {
	  ((ofstream*)_out)->close();
	  delete _out;
	  _out = NULL;
	  _firstTime=false;
	}
}

void myLog::printArgv(int loglvl, int argc, char *argv[]) {
  LOG(loglvl,<<"argv =");
  
  for (int i=0;i<argc;++i)
    LOG(loglvl,<<" \""<<argv[i]<<"\"");
  LOG(loglvl,<<endl);
  
}
