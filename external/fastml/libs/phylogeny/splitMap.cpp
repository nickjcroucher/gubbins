// $Id: splitMap.cpp 962 2006-11-07 15:13:34Z privmane $

#include "splitMap.h"
#include <algorithm>
using namespace std;

int splitMap::add(const split & in) { // add a split and return it's new count.
	return(_map[in]=_map[in]+1);
}

class valCmp {
public:
	bool operator()(const pair<split,int> & elem1, const pair<split,int> & elem2) {
		return (elem1.second > elem2.second);
	}
};

vector<pair<split,int> > splitMap::sortSplits() const{
	vector<pair<split,int> > svec(_map.size());
	partial_sort_copy(_map.begin(),_map.end(),svec.begin(),svec.end(),valCmp());
	return svec;
}

int splitMap::counts(const split& in) const {
	mapSplitInt::const_iterator i(_map.find(in));   
	if (i==_map.end()) return 0;
	return i->second;
}

void splitMap::print(ostream& sout) const {// default cout.
	for (mapSplitInt::const_iterator i = _map.begin(); i != _map.end();++i) {
      sout << i->second<<"\t"<<i->first;
    }
	sout <<endl;
}


ostream& operator<< (ostream &sout,  const splitMap& split_map) {
	split_map.print(sout);
	return sout;
}

/*splitMap::reverse_mapSplitInt splitMap::reverse() const
{
  reverse_sMap_t rmap;
	for (sMap_t::const_iterator i=_map.begin(); i!=_map.end();++i)
		rmap.insert(rMapPair_t(i->second,i->first));
	return rmap;
}
*/
