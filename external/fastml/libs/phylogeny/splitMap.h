// $Id: splitMap.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___SPLITMAP
#define ___SPLITMAP

#include "definitions.h"
#include "split.h"
#include <map>
using namespace std;

// splitMap is a map of split to integers used for counting the occurences of each split.
// Questions we want the class to be able to answer:
// 1. What is the occurence a specific split.
// 2. what is the most common split
// 3. Sort the splits according to their frequency.

class splitMap  {
// public:
// typedef pair<int,const split> rMapPair_t;
// typedef multimap<const int,const split> reverse_sMap_t;
// typedef multimap<int,split> reverse_sMap_t;
// reverse_sMap_t reverse() const ;
public:
  explicit splitMap(){};  // empty constractor
  int add(const split & in); // return the new frequency.
  int counts(const split& in) const;  // counts the number of occurances
  void print(ostream& sout = cout) const;
  vector<pair<split,int> > sortSplits() const;
private:
      
  typedef map<split,int> mapSplitInt;
  mapSplitInt _map;
};

ostream& operator<< (ostream &sout,  const splitMap& split_map);
#endif

