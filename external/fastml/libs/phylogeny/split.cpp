// $Id: split.cpp 962 2006-11-07 15:13:34Z privmane $

#include "split.h"
#include <cassert>
#include <algorithm>
using namespace std;

// there are always two options. Either the active set is _set[0] or _set[1].
// this depends on the parameter _reverse.
// The "1" will always be in the active set.
// so, for example consider the leaves [0,1,2] (_max = 3).
// The split {}{0,1,2} can be represented by both the empty split {} or the
// {0,1,2} split. Because the {0,1,2} split contains the "0" - this will be the active split.
// so we set _set[0] to be empty, and in _set[1] which is the active one (_reverse = true)
// we insert the leaves.
split::split (const int max): _max(max), _reverse(true){
	for(int j=0;j<max;++j) {
		_set[1].insert(j);
	}
}

// isMember searches for the key in the active set.
bool split::isMember(const int key) const {
	return(_set[_reverse].find(key)!=_set[_reverse].end());
}


void split::reverseMembership(const int key){
  assert(key<_max && key >= 0);

  // where is the key now
  // if the key is member, than in = _reverese;
  // Otherwise in = !_reverse
  bool in =(isMember(key))?_reverse:!_reverse; 

 _set[in].erase(key);
 _set[!in].insert(key);
  if (key==0)		// if we add "0", we need to reverse the split
    reverse();
};


int split::size() const  {
	int tmp = _set[_reverse].size();
	return (tmp<_max-tmp?tmp:_max-tmp);
}

void split::print(ostream& sout) const{ //  = cout
	sout <<"size ="<<size()<<"   ";
	set<int>::const_iterator i;
    for (i=_set[_reverse].begin();i != _set[_reverse].end();++i)
		sout << *i << " ";
    sout <<" | ";
    for (i=_set[!_reverse].begin();i != _set[!_reverse].end();++i)
		sout << *i << " ";
    sout << endl;
}

bool split::lessThen(const split& other) const{
    return(_set[_reverse]<other._set[other._reverse]);
}

bool split::compatible(const split& other) const {
    set<int>::const_iterator i     (_set[_reverse].begin());
    set<int>::const_iterator i_end (_set[_reverse].end());
    set<int>::const_iterator j     (other._set[other._reverse].begin());
    set<int>::const_iterator j_end (other._set[other._reverse].end());
	return (includes(i,i_end,j,j_end) || includes(j,j_end,i,i_end));
}

void split::reverse(){		// actualy reverse membership in the set
	_reverse=!_reverse;
 }

bool operator<(const split& a, const split& b) {
	return(a.lessThen(b));
}

ostream& operator<< (ostream &sout,  const split& split) {
	split.print(sout);
	return sout;
}


