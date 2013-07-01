// $Id: phylipFormat.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___PHYLIP_INTERLEAVED_FORMAT
#define ___PHYLIP_INTERLEAVED_FORMAT

#include "definitions.h"
#include "sequenceContainer.h"

class phylipSequentialFormat {
public:
	static sequenceContainer read(istream &infile, const alphabet* alph);
	static void write(ostream &out, const sequenceContainer& sd,
		const int numOfPositionInLine = 50,
		const int spaceEvery = 10);
	//readUnAligned: the input sequences do not need to be aligned (not all sequences are the same length).
	static sequenceContainer readUnAligned(istream &infile, const alphabet* alph);
};

#endif

/* EXAMPLE OF PHYLIP FORMAT (sequential):

6   128
Langur     KIFERCELAR TLKKLGLDGY KGVSLANWVC LAKWESGYNT EATNYNPGDE
		   STDYGIFQIN SRYWCNNGKP GAVDACHISC SALLQNNIAD AVACAKRVVS
		   DQGIRAWVAW RNHCQNKDVS QYVKGCGV
Baboon     KIFERCELAR TLKRLGLDGY RGISLANWVC LAKWESDYNT QATNYNPGDQ
           STDYGIFQIN SHYWCNDGKP GAVNACHISC NALLQDNITD AVACAKRVVS
		   DQGIRAWVAW RNHCQNRDVS QYVQGCGV
Human      KVFERCELAR TLKRLGMDGY RGISLANWMC LAKWESGYNT RATNYNAGDR
           STDYGIFQIN SRYWCNDGKP GAVNACHLSC SALLQDNIAD AVACAKRVVR
           DQGIRAWVAW RNRCQNRDVR QYVQGCGV

*/

