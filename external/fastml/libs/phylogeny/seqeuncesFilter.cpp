#include "seqeuncesFilter.h"
#include "nucleotide.h"

seqeuncesFilter::~seqeuncesFilter()
{}

void seqeuncesFilter::removeSequencesWithStop(sequenceContainer & sc, codon & alpha)
{

	//going over al seqeunces
	for (int i = 0; i < sc.numberOfSeqs();++i) {
		int id = sc.placeToId(i);
		//going over all sequence len
		for (int j = 0; j < sc.seqLen();++j) { 	
			//remove seqeunces with stop data not in the middle
			if ((j != sc.seqLen()-1) && (alpha.isStopCodon(sc[id][j])))
			{ 
				LOG(4, <<"removing sequence = "<<sc.name(id)<<" : STOP codon in the middle of the reading frame!"<<endl);
				sc.remove(id);
				i--;
				break;
			}
		}
	}
}

void seqeuncesFilter::removeSequencesWithMissingData(sequenceContainer & sc)
{

	//going over al seqeunces
	for (int i = 0; i < sc.numberOfSeqs(); ++i) 
	{	
		//going over all sequence len
		for (int j = 0; j < sc.seqLen(); ++j) 
		{ 
			int id = sc.placeToId(i);
			//remove seqeunces with unkonwn data
			if (sc[id][j] == sc.getAlphabet()->unknown())
			{ 
				sc.remove(id); 
				i--;
				break;
			}
		}
	}
}

void seqeuncesFilter::removeSequencesWithMissingDataAndStop(sequenceContainer & sc, codon & alpha)
{

	//going over al seqeunces
	for (int i = 0; i < sc.numberOfSeqs(); ++i) {
		int id = sc.placeToId(i);
		//going over all sequence len
		for (int j = 0; j < sc.seqLen();++j) { 	
			//remove seqeunces with stop data not in the middle or missing data
			if ((j != sc.seqLen()-1) && (sc[id][j] == sc.getAlphabet()->unknown() || alpha.isStopCodon(sc[id][j])))
			{	
				
				sc.remove(id);
				i--;
				break;
			}
		}
	}

}


void seqeuncesFilter::removeSequencesNotStartWithATG(sequenceContainer & sc, codon & alpha)
{
	amino aa;
	//going over al seqeunces
	for (int i = 0; i < sc.numberOfSeqs();++i) {
		int id = sc.placeToId(i);
		int in_first = codonUtility::aaOf(sc[id][0], alpha);
		if (in_first != aa.fromChar('M'))
		{
				LOG(4, <<"removing sequence = "<<sc.name(id)<<" : not starting with ATG!"<<endl);
				sc.remove(id);
				i--;
		}
	}
}

void seqeuncesFilter::removeSequencesNotStartWithInitiationCodons(sequenceContainer & sc,codon & alpha)
{
	for (int i = 0; i < sc.numberOfSeqs();++i) {
		int id = sc.placeToId(i);
		int in_first = sc[id][0];
		if(!alpha.isInitiationCodon(in_first)){
			LOG(4, <<"removing sequence = "<<sc.name(id)<<" : not starting with initiation codon!"<<endl);
			sc.remove(id);
			i--;
		}
	}
}


void seqeuncesFilter::removeSequencesWithGapsAccordingRef(sequenceContainer & sc,int precent,string refName)
{
	int refID = sc.getId(refName);
	Vint seqToRemove;
	//going over all position in reference seqeunce
	for (int pos = 0; pos < sc[refID].seqLen(); pos++)
	{	
		
		//check if the pos is gap
		if (sc[refID][pos] == sc.getAlphabet()->gap())
			//going over all other seqeunces to compute the precents of gaps
		{
			cout<<pos<<" ";
			seqToRemove.clear();
			MDOUBLE numOfSeqWithOutGap = 0;
			cout<<sc.numberOfSeqs()<<" ";
			for (int i = 0; i < sc.numberOfSeqs(); i++)
			{
				
				int id = sc.placeToId(i);
				if  (sc[id][pos] != sc.getAlphabet()->gap())
				{
					numOfSeqWithOutGap++;
					seqToRemove.push_back(id);
				}
			}
			cout<<seqToRemove.size()<<endl;
			if ((100 * ((sc.numberOfSeqs() - numOfSeqWithOutGap)/sc.numberOfSeqs())) > precent)
			{	
				for (int j = 0; j < seqToRemove.size(); j++){
					sc.remove(seqToRemove[j]);
				}
			
			}
		}
	}
}

//removes all sequences that are shorter than lowerBound and longer than upperBound
void seqeuncesFilter::removeShortAndLongSequences(sequenceContainer & sc, int lowerBound, int upperBound)
{
	const alphabet* pAlph = sc.getAlphabet();
	//going over al seqeunces
	for (int seq = 0; seq < sc.numberOfSeqs(); ++seq) 
	{
		int id = sc.placeToId(seq);
		//checking sequence length
		int seqLen = sc[id].seqLenSpecific();
		if ((seqLen < lowerBound) || (seqLen > upperBound))
		{
			cerr<<"removing sequence: "<<sc.name(id)<<" sequence Length = "<<seqLen<<endl;
			sc.remove(id);
			--seq;
		}
	}
}

//removes all sequences that have inserts in which most other sequences (> percent) have gaps.
//in case refName is given: check only positions in which the reference sequence has gaps.
//The remained sequences are stored in newSc.
void seqeuncesFilter::removeSequencesWithInserts(sequenceContainer & newSc,const sequenceContainer & sc,int percent, const string& refName, string outFileName)
{
	if (outFileName.empty())
		outFileName = "removedSequences" + double2string(percent) + ".txt";
	ofstream outF(outFileName.c_str());
	int refID;
	if (!refName.empty())
		refID = sc.getId(refName);
	Vint seqToAdd(sc.numberOfSeqs(), 1);//1== add the sequence to newSc. 0 = don't add.
	//going over all position (in reference seqeunce if given)
	for (int pos = 0; pos < sc.seqLen(); ++pos)
	{	
		
		if (!refName.empty())
		{	//don't remove this position if it isn't gap in the refSeqeunce
			if (sc[refID][pos] != sc.getAlphabet()->gap())
				continue;
		}
		Vint seqToRemove; //holds the ids of sequences without gaps in the current positions
		//going over all seqeunces to compute the percent of gaps
		MDOUBLE numOfSeqWithGap = 0;
		for (int i = 0; i < sc.numberOfSeqs(); i++)
		{
			int id = sc.placeToId(i);
			if  (sc[id][pos] != sc.getAlphabet()->gap())
			{
				seqToRemove.push_back(id);
			}
			else
				numOfSeqWithGap++;
		}
		//outF<<"POS "<<pos<<" seqWithGaps = "<<numOfSeqWithGap<<" seqWithoutGaps = "<<sc.numberOfSeqs() - numOfSeqWithGap<<endl;
		//in case most sequences have gaps in that position: remove the sequences that have inserts at that position
		MDOUBLE percentGapsinPos = 100.0 * (numOfSeqWithGap / sc.numberOfSeqs());
		if (percentGapsinPos > percent)
		{
            //outF<<"removing sequences: ";
			for (int j = 0; j < seqToRemove.size(); j++)
			{
				int x = seqToRemove[j];
				seqToAdd[seqToRemove[j]] = 0;
				outF<<sc.name(sc.placeToId(x))<<endl;
			}
			outF<<endl;
		}
	}
	
	
	for (int i=0; i<seqToAdd.size(); i++)
	{	
		if (seqToAdd[i] == 1)
		{
			int id = sc.placeToId(i);
			newSc.add(sc[id]);
		}
	}
	outF.close();
}

void seqeuncesFilter::removeSequencesNotDivisableBy3(sequenceContainer & sc)
{
	nucleotide nucAlph;
	for (int i = 0; i < sc.numberOfSeqs();++i) 
	{
		int id = sc.placeToId(i);
		int seqL = sc[id].seqLen();
		if ((seqL % 3) != 0)
		{
				LOG(4, <<"removing sequence = "<<sc.name(id)<<" : nucleotide sequence length is not divisable by 3!"<<endl);
				sc.remove(id);
				--i;
		}
	}
}
