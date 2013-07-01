// $Id: codon.h 5975 2009-03-17 08:00:37Z rubi $
#ifndef ____CODON
#define ____CODON

#include <cassert>
#include "definitions.h"
#include "errorMsg.h"
#include "someUtil.h"
#include "alphabet.h"
#include "geneticCodeHolder.h"
#include <map>
class codon;

class sequenceContainer;
class codonUtility {
public:
	enum diffType {equal =0, tr, tv, twoTrs, twoTvs ,trtv, threesub};
	static diffType codonDiff(const int c1, const int c2, const codon &cod);
	static diffType codonDiff(const int c1, const int c2) {return _trtvDiff[c1][c2];}

	enum replacementType {sameCodon=0, synonymous, non_synonymous};
	static replacementType codonReplacement(const int c1, const int c2, const codon &cod);
	static replacementType codonReplacement(const int c1, const int c2) {return _synNonsynDiff[c1][c2];}

	enum nucDiffPlaceType {A1=0, A2, A3,C1, C2, C3, G1,G2,G3,T1,T2,T3, EQUAL, MUL_SUB};
	static nucDiffPlaceType nucDiffPlace(const int fromCodon, const int targetCodon, const codon &cod);
	static nucDiffPlaceType nucDiffPlace(const int fromCodon, const int targetCodon) {return _nucDiffPlace[fromCodon][targetCodon];}

	enum nucsDiffType {AC=0, AG, AT, CG, CT, GT, SAME, DIFF}; //The difference between two codons: For exampe nucsDiff(ACT, ACG) returns GT. DIFF = more than one change. 
	static nucsDiffType nucsDiff(const int fromCodon, const int targetCodon, const codon &cod);
	static nucsDiffType nucsDiff(const int fromCodon, const int targetCodon) {return _nucsDiff[fromCodon][targetCodon];}

	static int aaOf(const int c1, const codon &cod);
	static void initSubMatrices(const codon& cod);
	
	//returns the number (codonCounter) and frequency (codonUsage) of each codon in the sequnece container
	static void getCodonUsage(const sequenceContainer& sc, Vint& codonCounter, Vdouble& codonUsage);
	static void readCodonUsage(const string& codonUsageFileName, Vdouble& codonUsage,const codon &inCodonAlpa);
	//calculates the CAI for the whole MSA and for each position. 
	//The calculation is based on a pre-calculated codonUsage vector.
	static MDOUBLE calcCodonAdaptationIndex(const sequenceContainer& sc, const Vdouble& codonUsage, Vdouble& cai4site);

private:
	static vector<vector<diffType> > _trtvDiff;
	static vector<vector<replacementType> > _synNonsynDiff;
	static vector<vector<nucDiffPlaceType> > _nucDiffPlace;
	static vector<vector<nucsDiffType> > _nucsDiff;
};


class codon : public alphabet {
public:
	explicit codon(); //default constructor: reads "nuclearCode.txt"
	explicit codon(const geneticCodeString& matrixFileString);
	virtual ~codon() {}
  //	explicit codon( codon& other);
	codon& operator=(const codon& other);
	virtual alphabet* clone() const { return new codon(*this); }
	void readMatrixFromFile(const string& matrixFileName);
	const map <string,string> & geneticCode()const {return _geneticCode;}
	int unknown() const  {return 64;}
	int gap() const  {return -1;}
	int size() const {return _alphabetSize;} // 3 stop codon excluded
	int stringSize() const {return 3;} // 3 letter code.
	vector<int> fromString(const string& str) const;
	bool isStopCodon(const int in_id) const; 
	bool isStopCodon(const string& str) const {return isStopCodon(fromChar(str));}; 
	bool isInitiationCodon(const int in_id) const; 
	bool isInitiationCodon(const string& str) const {return isInitiationCodon(fromChar(str));}; 
	int fromChar(const string& s, const int pos=0) const;
	string fromInt(const int in_id) const;
	// "specific" here is not unknown, nor ambiguity, nor gap (for example, for nucleotides it will true for A,C,G, or T).
	bool isSpecific(const int id) const {return (id>=0 && id < size());}


	
  int relations(const int charInSeq, const int charToCheck) const{
		if (charInSeq == -1) {
			errorMsg::reportError("gaps in the sequences. Either change gaps to ? or remove gap positions");
		}
		else if (charInSeq == unknown()) return 1;
		else if (charInSeq == charToCheck) return 1;		
		if (charInSeq >= _alphabetSize)
		{
			string err= "";
			err+="charInSeq = ";
			err += int2string(charInSeq);
			err+= " _alphabetSize = ";
			err+=int2string(_alphabetSize);
			errorMsg::reportError(err);
		}
		assert(charInSeq < _alphabetSize);
		return 0;
	}
private:
  void init(const geneticCodeString& matrixFileString);
private:
	map <string,string> _geneticCode; //key - codon, value - amino acid
	map <string,int> _codon2Int;//key string of codon int= integer value of codon
	map <int,string> _initiationIndex2codon;//key: integer value of codon; value: string of initiation codon. the keys is an integer so that the value of the init codon can be found  
	int _alphabetSize;
};




#endif
