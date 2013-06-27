// $Id: readTree.cpp 5525 2008-12-19 20:17:05Z itaymay $

#include "definitions.h"
#include "errorMsg.h"
#include "someUtil.h"
#include "readTree.h"
#include <iostream>
using namespace std;





// forward declarations

//----------------------------------------------------------------------------------------------
//	about reading tree topology from files:
//	usually a tree topology is represented by a line like this
//	(((Langur:0.8,Baboon:0.55):0.3,Human:0.44):0.5,Rat:0.02,(Cow:0.2,Horse:0.04):0.03);
//	the syntax of such a line is (part, part, part, part)
//	where part is either (part,part, part, ...):distace or name:distance
// or without the distance!
//	it should notice that the tree is unrooted.
//	if we look at the above file format, one can notice that the number of comas (",") is 
//	always one less than the number of leaves (synonyms for leaves are OTUs and external nodes)
//	the function GetNumberOfLeaves counts the numnber of comas and returns the number of leaves.
//	in the example below there are 6 leaves.

//*******************************************************************************
// constructors
//*******************************************************************************





vector<char> PutTreeFileIntoVector(istream &in) {
	vector<char> tree_contents;
	bool endWithDotComa = false;
	char chTemp;
	while (( !in.eof()) && (tree_contents.size() < MAX_FILE_SIZE))
	{
		in.get(chTemp);
#ifdef WIN32
		if (chTemp == -52) return tree_contents; //tal addition.
#endif
		if ( !isspace( chTemp ) )
			tree_contents.push_back(chTemp);
		if (chTemp == ';') {
			endWithDotComa = true;
			break;
		}
	}

	if (tree_contents.size() >= MAX_FILE_SIZE) {
		vector<string> err;
		err.push_back("Error reading tree file. The tree file is too large");
		errorMsg::reportError(err,1); // also quit the program
	}
	if (endWithDotComa == false) tree_contents.clear(); // remove junk from the last ; till the end of the file.
	return tree_contents;
}




int GetNumberOfLeaves(const vector<char> &tree_contents) {
	int iCommasCounter = 0;
	vector<char>::const_iterator itCurrent = tree_contents.begin();
	for ( ; itCurrent != tree_contents.end(); ++itCurrent ) {
		if (*itCurrent==COMMA)
			++iCommasCounter;
	}
	return ++iCommasCounter; //#leaves is always one more than number of comas
}

int GetNumberOfInternalNodes(const vector<char> &tree_contents) {
	int iCloseCounter = 0;
	vector<char>::const_iterator itCurrent = tree_contents.begin();
	for ( ; itCurrent != tree_contents.end(); ++itCurrent ) {
		if (*itCurrent==CLOSING_BRACE) ++iCloseCounter;
		if (*itCurrent==CLOSING_BRACE2) ++iCloseCounter;
	}
	return iCloseCounter; //number of HTUs is always the number of ")"
}


bool verifyChar(vector<char>::const_iterator &p_itCurrent, const char p_cCharToFind) {
	if ( (*p_itCurrent)==p_cCharToFind ) 	return true;
	return false;
}




// IsAtomicPart decides whether we will now read a taxa name (return true),
// or read an OPENING_BRACE which will say us, that we will read a complicated strucure.
bool IsAtomicPart(const vector<char>::const_iterator p_itCurrent) {
	if ( (*p_itCurrent)==OPENING_BRACE ) return false;
	else if ( (*p_itCurrent)==OPENING_BRACE2 ) return false;
	return true;
}

//-----------------------------------------------------------------------------
// there are 2 options for the tree format.
// either (name1:0.43, name2: 0.45 , (name3 : 2 , name 4: 5) : 3.332)
// or without the distances (name1, name2 , (name3 , name4) )
// here we return true if the tree file is with the distance, or false, if the tree file
// has not distances. 
// if distances exist: after the name there will always be a colon
// if distance exist, also move the iterator, to the beggining of the number
//-----------------------------------------------------------------------------
bool DistanceExists(vector<char>::const_iterator& p_itCurrent) {

	if ((*p_itCurrent)==COLON ) 	{
		++p_itCurrent;
		return true;
	}
	return false;
}

void clearPosibleComment(vector<char>::const_iterator& p_itCurrent) {
  if ((*p_itCurrent)=='[' ) {
	while (*(++p_itCurrent) != ']');
	++p_itCurrent;				// move over "]"
  }
}

string readPosibleComment(vector<char>::const_iterator& p_itCurrent) {
  string comment = "";

  if ((*p_itCurrent)=='[' ) 
  {
	vector<char>::const_iterator tmp= (p_itCurrent+1);
	if ((*tmp++)=='&' &&
		(*tmp++)=='&' &&
		(*tmp++)=='N' &&
		(*tmp++)=='H' &&
		(*tmp++)=='X') // see http://www.genetics.wustl.edu/eddy/forester/NHX.pdf
	// [&&NHX...]
	{			
		p_itCurrent += 5;
		while (*(++p_itCurrent) != ']')
		{
		  comment += *(p_itCurrent);
		}
		++p_itCurrent;				// move over "]"
	}
	else // [...]
	{
		// Skip over the text in []
		++p_itCurrent;
		while (*(p_itCurrent) != ']')
            ++p_itCurrent;				
		++p_itCurrent;	// move over "]"
		
	}
  }
  if (comment.size())
	LOG(10,<<"comment ="<<comment<<endl);

  return comment;
}



MDOUBLE getDistance(vector<char>::const_iterator &p_itCurrent) {
	string sTempNumber;
	for ( ; isdigit(*p_itCurrent) || (*p_itCurrent)==PERIOD || (*p_itCurrent)=='E'|| (*p_itCurrent)=='e'|| (*p_itCurrent)=='-' || (*p_itCurrent)=='+'; ++p_itCurrent)	
		sTempNumber += (*p_itCurrent);
	MDOUBLE dDistance = string2double(sTempNumber);
	return dDistance;
}





