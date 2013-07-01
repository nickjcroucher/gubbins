// $Id: distribution.cpp 2709 2007-11-19 14:49:21Z itaymay $

#include "distribution.h"
#include "errorMsg.h"

distribution::~distribution(){}
// this must be here. see Effective c++ page 63 (item 14, constructors, destructors, 
// assignment

void distribution::change_number_of_categories(int in_number_of_categories)
{
	errorMsg::reportError("not implemented: distribution::change_number_of_categories()!");
}