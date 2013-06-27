#include "samplingSequences.h"
#include "logFile.h"
#include "talRandom.h"


sampleSequences::sampleSequences(sequenceContainer &sc){
	_sc = sc;
}

sequenceContainer sampleSequences::removeSequences(sequenceContainer &sc){
	int noOfSeq = sc.numberOfSeqs();
	int gap = sc.getAlphabet()->gap();
	int unknown = sc.getAlphabet()->unknown();
	bool seqToAdd;
	int n =0;
	sequenceContainer newSc;
	for (int i=0;i<noOfSeq;i++){
		seqToAdd = true;
		for (int j=0;j<sc[i].seqLen();j++){
			if ((sc[i][j]== gap) || (sc[i][j]== unknown ) || (sc[i].seqLen() != 297)){
				seqToAdd = false;
			}
		}
		if (seqToAdd == true) {
			sequence add = sc[i];
			sequence sc(add);
			sc.setID(n);
			n++;
			newSc.add(sc);
		}
	}
	return newSc;
}


void sampleSequences::printDistances(){
	for (int i=0;i< _distances.size();i++){
		for (int j=0;j<_distances[i].size();j++){
			cout<<_distances[i][j]<<" ";
		}
		cout<<endl;
	}
}

void sampleSequences::setDistance(int i,int j,MDOUBLE dist){
	(i<j ? _distances[i][j-i] :_distances[j][i-j]) = dist;
}

MDOUBLE sampleSequences::getDistance(int i,int j){	
	return (i<j ? _distances[i][j-i] :_distances[j][i-j]);
}


sequenceContainer sampleSequences::sampleFarthestSequences(int n, distanceMethod *dm){
	_sc.removeIdenticalSequences();
	if (n >= _sc.numberOfSeqs()){
		cerr<<"Number of sequences to sample is bigger than the origin number of sequences so the all sequences were chosen in sampleSequences::sampleFarthestSequences"<<endl;
		return _sc;
	}

	int numberOfSeq = _sc.numberOfSeqs();
	_distances.resize(numberOfSeq);
	int i;
	for (i=0;i<numberOfSeq;i++)
		_distances[i].resize(numberOfSeq-i);

	for (i=0;i<numberOfSeq;i++){
		for(int j=i;j<numberOfSeq;j++){	
			int id1 = _sc.placeToId(i);
			int id2 = _sc.placeToId(j);

			setDistance(i,j,dm->giveDistance(_sc[id1],_sc[id2],NULL));
		}	
	}

	sequenceContainer newSc;
	vector<int> sampled;
	sampled.push_back(0);//to change
	int id = 0;
	int p = _sc.placeToId(0);
	sequence sc(_sc[p]);
	sc.setID(id++);
	newSc.add(sc);
	while (newSc.numberOfSeqs()<n){
		int i = findNextSeq(sampled);
		p = _sc.placeToId(i);		
		sequence sc(_sc[p]);
		sc.setID(id);
		newSc.add(sc);
		id++;		
		sampled.push_back(i);
	}
	return newSc;
}

int sampleSequences::findNextSeq(vector<int> &sampled){	
	MDOUBLE max = 0,min;
	int seqi = -1;
	for(int i=0;i< _sc.numberOfSeqs();i++){
		min=10000;//to update
		for (int j=0;j<sampled.size();j++){
			if (getDistance(i,sampled[j])<min)
				min = getDistance(i,sampled[j]);
		}
		if (max<min){
			max=min;
			seqi = i;
		}		
	}
	
	if (seqi>_sc.numberOfSeqs() ||seqi<0){
		errorMsg::reportError("Error in sampleSequences::findNextSeq");
	}
	return seqi;
}

//sequenceContainer sampleSequences::sampleRandomSequences(int seqNum)
//{
//	if (seqNum > _sc.numberOfSeqs())
//		errorMsg::reportError("sampleSequences::sampleRandomSequences(): the number of requested seqeuences is larger than the number of sequences in the MSA");
//	sequenceContainer newSc(_sc);
//	while (newSc.numberOfSeqs() > seqNum)
//	{
//		int seqPlaceToRemove = talRandom::giveIntRandomNumberBetweenZeroAndEntry(newSc.numberOfSeqs());
//        newSc.remove(newSc.placeToId(seqPlaceToRemove));
//	}
//	return newSc;
//}


sequenceContainer sampleSequences::sampleRandomSequences(int seqNum)
{
	if (seqNum > _sc.numberOfSeqs())
		errorMsg::reportError("sampleSequences::sampleRandomSequences(): the number of requested seqeuences is larger than the number of sequences in the MSA");
	sequenceContainer newSc;
	Vint vec2Add(_sc.numberOfSeqs(),0);
	int n = 0;
	while (n < seqNum)
	{
		int seqPlaceToAdd = talRandom::giveIntRandomNumberBetweenZeroAndEntry(_sc.numberOfSeqs());
		if (vec2Add[seqPlaceToAdd] == 0){
			vec2Add[seqPlaceToAdd] = 1;
			n++;
		}
        
	}
	for (int i = 0; i<vec2Add.size();i++){
		if (vec2Add[i] == 1)
			newSc.add(_sc[i]);
	}
	return newSc;
}
//sequenceContainer sampleSequences::sampleRandomCharacters(int seqLen)
//{
//	if (seqLen > _sc.seqLen())
//		errorMsg::reportError("sampleSequences::sampleRandomCharacters(): the requested sequence length is larger than the number of characters in the MSA");
//	Vint posToRemove(_sc.seqLen(),1);
//	//first create a vector with seqLen positions to be sampled in the begining of the vector
//	for (int i = 0; i < seqLen; ++i)
//		posToRemove[i] = 0;
//	//then randomly swap the positions in posToRemove.
//	//The end result is a random vector with the positions to remove marked with '1'
//	int swapNum = _sc.seqLen() * 10;
//	for (int x = 0; x < swapNum; ++x)
//	{
//		int pos1 = talRandom::giveIntRandomNumberBetweenZeroAndEntry(_sc.seqLen());
//		int pos2 = talRandom::giveIntRandomNumberBetweenZeroAndEntry(_sc.seqLen());
//		int tmp = posToRemove[pos1];
//		posToRemove[pos1] = posToRemove[pos2];
//		posToRemove[pos2] = tmp;
//	}
//	
//	sequenceContainer newSc(_sc);
//	newSc.removePositions(posToRemove);
//	return newSc;
//}


sequenceContainer sampleSequences::sampleRandomCharacters(int seqLen)
{
	if (seqLen > _sc.seqLen())
		errorMsg::reportError("sampleSequences::sampleRandomCharacters(): the requested sequence length is larger than the number of characters in the MSA");
	sequenceContainer newSc(_sc);

	while (newSc.seqLen() > seqLen)
	{
		Vint posToRemove(newSc.seqLen(),0);
		int seqPlaceToRemove = talRandom::giveIntRandomNumberBetweenZeroAndEntry(newSc.seqLen());
		posToRemove[seqPlaceToRemove] = 1;
		newSc.removePositions(posToRemove);        
	}
	return newSc;
}
