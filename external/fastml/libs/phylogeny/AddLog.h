// $Id: AddLog.h 962 2006-11-07 15:13:34Z privmane $

// version 1.00
// last modified 2 Nov 2002

#ifndef __AddLog_h
#define __AddLog_h

#include <iostream>
using namespace std;

class tAddLog_Precompute {
  public:

  tAddLog_Precompute();
  ~tAddLog_Precompute();

  double AddLog( double x, double y );
  
private:
  static const int D_LOGADD; // = 50;   // y/x < 1e-D discard
  static const int G_LOGADD;// = 500;  // step function look-up every 1/G
  static int d_logadd;

  double *logaddf;
};

extern tAddLog_Precompute AddLogData;

inline
double
AddLog(double x, double y ){
  return AddLogData.AddLog(x, y);
}

inline double
tAddLog_Precompute::AddLog(double x, double y ){
  if (x < y)  {
    double dummy = x;
    x = y;
    y = dummy;
  }

#ifdef notdef  
  return x + log(1 + exp(y-x));
#endif
  
  double z = (x-y)*G_LOGADD;
  int i = int(z);
  if( i < d_logadd ) x += ((i+1-z)*logaddf[i] + (z-i)*logaddf[i+1]);
  return x;
}

#endif


/*
Folks,

In many of our program we use the AddLog procedure that compute the sum of
two numbers in log form. Gill spent some time investigating faster versions
of this procedure, which gave him 3-4 fold speedup on his program. Attached
is my re-packaging of his solution. I think it will be useful in some of the
code we use.

-Nir
*/
