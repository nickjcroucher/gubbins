// $Id: uniDistribution.cpp 2711 2007-11-19 14:49:54Z itaymay $

#include "uniDistribution.h"
#include "errorMsg.h"


void uniDistribution::change_number_of_categories(int in_number_of_categories)
{
	if (in_number_of_categories != 1)
		errorMsg::reportError("error in uniDistribution::change_number_of_categories() - number of categories is not 1");
}
