// $Id: logFile.h 6067 2009-04-14 19:12:28Z itaymay $

#ifndef ___LOG
#define ___LOG


#include <string>
#include <iostream>
#include <fstream>

using namespace std;
			
class myLog {
public:
	static int LogLevel() { return _loglvl;}
	static ostream& LogFile(void) {
		if (_out == NULL) return cerr;
		return *_out;
	}

	static void setLogLvl(const int newLogLvl) {_loglvl = newLogLvl;}
	static void setLogOstream(ostream* out) {_out = out;}

	// this function is problematic, because it issue a call to NEW
	// which because the function is static - cannot be deleted.
	// but, this will not effect the program, because there is only
	// 1 instance of _out and it will be released anyway in the end of the program.
	static void setLog(const string logfilename, const int loglvl);
	static void endLog(void);
	static void printArgv(int loglvl, int argc, char *argv[]) ;
private:
	static ostream* _out;
	static int _loglvl;
	static bool _firstTime;
};

#ifdef LOG
#undef LOG
#endif
		

#define LOG(Lev, ex) { if( Lev <= myLog::LogLevel() ) myLog::LogFile() ex; }
#define LOGnOUT(Lev, ex) { if( Lev <= myLog::LogLevel() ) {myLog::LogFile() ex; cerr ex; }}
#define LOGDO(Lev, ex) { if( Lev <= myLog::LogLevel() ) ex; }


#endif



