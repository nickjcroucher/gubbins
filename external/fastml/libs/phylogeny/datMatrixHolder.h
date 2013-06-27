// $Id: datMatrixHolder.h 5804 2009-01-20 09:18:05Z adido $

#ifndef ___DATMATRIXHOLDER
#define ___DATMATRIXHOLDER

#include <string>
using namespace std;

// THIS CONSTRUCT IS USED TO KEEP A STRING THAT IS THE AA SUBSTITUTION MATRIX
// THE datMatrixString IS TO BE USED WHENEVER WE USE ONE OF THE BUILD-IN AA SUBSTITUTION MATRICES.

class datMatrixString {
public:
  const string Val;
  explicit datMatrixString(const char * str): Val(str){};
};

class datMatrixHolder {
public:
  static const datMatrixString cpREV45;
  static const datMatrixString dayhoff;
  static const datMatrixString jones;	// This is JTT
  static const datMatrixString mtREV24;
  static const datMatrixString wag;
  static const datMatrixString HIVb;
  static const datMatrixString HIVw;
  static const datMatrixString lg;
  static const datMatrixString empiriCodon; //This is the empirical matrix for codon by gina and adrian
};

#endif	// ___DATMATRIXHOLDER
