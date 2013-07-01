#ifndef SAMPLE_SEQUENCES_H
#define SAMPLE_SEQUENCES_H

#include "definitions.h"
#include "distanceMethod.h"
#include "sequenceContainer.h"
#include "pDistance.h"


class sampleSequences{
public:
	explicit sampleSequences(sequenceContainer &sc);	
	virtual ~sampleSequences() {};

	sequenceContainer sampleFarthestSequences(int n, distanceMethod *dm);
	//sampleRandomSequences: samples seqNum sequences from the sequence container
	sequenceContainer sampleRandomSequences(int seqNum);
	//sampleRandomCharacters: samples seqLen characters from the sequenceContainer
	sequenceContainer sampleRandomCharacters(int seqLen);


private:
	int findNextSeq(vector<int> &sampled);
	void setDistance(int i,int j,MDOUBLE dist);
	MDOUBLE getDistance(int i,int j);
	void removeSequenceWithGap();
	sequenceContainer removeSequences(sequenceContainer &sc);
	void printDistances();
private:
	VVdouble _distances;
	sequenceContainer _sc;
};
#endif
