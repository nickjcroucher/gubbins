// 	$Id: geneticCodeHolder.h 962 2006-11-07 15:13:34Z privmane $	

#ifndef ___GENMATRIXHOLDER
#define ___GENMATRIXHOLDER

#include <string>
using namespace std;

// THIS CONSTRUCT IS USED TO KEEP A STRING THAT IS THE AA SUBSTITUTION MATRIX
// THE datMatrixString IS TO BE USED WHENEVER WE USE ONE OF THE BUILD-IN AA SUBSTITUTION MATRICES.

class geneticCodeString {
public:
  const string Val;
  explicit geneticCodeString(const char * str): Val(str){};
};

class geneticCodeHolder {
public:
  static const geneticCodeString nuclearStandard;
  static const geneticCodeString nuclearEuplotid;
  static const geneticCodeString nuclearCiliate;
  static const geneticCodeString nuclearBlepharisma;
  static const geneticCodeString mitochondriaYeast;
  static const geneticCodeString mitochondriaVertebrate;
  static const geneticCodeString mitochondriaProtozoan;
  static const geneticCodeString mitochondriaInvertebrate;
  static const geneticCodeString mitochondriaFlatworm;
  static const geneticCodeString mitochondriaEchinoderm;
  static const geneticCodeString mitochondriaAscidian;
};

#endif	// ___GENMATRIXHOLDER
