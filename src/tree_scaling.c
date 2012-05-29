/*
 *  Wellcome Trust Sanger Institute
 *  Copyright (C) 2012  Wellcome Trust Sanger Institute
 *  
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "seqUtil.h"
#include "Newickform.h"
#include "tree_scaling.h"

void scale_branch_distances(newick_node * root_node, int number_of_filtered_snps)
{

	if (root_node->childNum == 0)
	{
		root_node->dist = (root_node->dist * number_of_filtered_snps);
	}
	else
	{
		root_node->dist = root_node->dist * number_of_filtered_snps;
		newick_child *child;
		child = root_node->child;
		while (child != NULL)
		{
			scale_branch_distances(child->node, number_of_filtered_snps);
			child = child->next;
		}
	}
}
