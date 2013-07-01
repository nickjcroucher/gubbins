#ifndef __SEQEUNCES_FILTER
#define __SEQEUNCES_FILTER


#include "definitions.h"
#include "sequenceContainer.h"
#include "codon.h"
#include "amino.h"
#include <string>
#include <fstream>
#include "fastaFormat.h"


using namespace std;

class seqeuncesFilter{

public:
	static void removeSequencesWithStop(sequenceContainer & sc,codon & alpha);
	static void removeSequencesWithMissingData(sequenceContainer & sc);
	//applied only to coding nucleotide seqeunces: remove sequence that are not divisable by 3.
	static void removeSequencesNotDivisableBy3(sequenceContainer & sc);
	static void removeSequencesWithMissingDataAndStop(sequenceContainer & sc,codon & alpha);
	static void removeSequencesNotStartWithATG(sequenceContainer & sc,codon & alpha);
	static void removeSequencesNotStartWithInitiationCodons(sequenceContainer & sc,codon & alpha);
	static void removeSequencesWithGapsAccordingRef(sequenceContainer & sc,int precent, string refName);
	static void removeSequencesWithInserts(sequenceContainer & newSc, const sequenceContainer & sc, int percent, const string& refName = "", string outFileName = "");


	//removes all sequences that are shorter than lowerBound and longer than upperBound
	static void removeShortAndLongSequences(sequenceContainer & sc, int lowerBound, int upperBound);
	virtual ~seqeuncesFilter();

};
#endif
