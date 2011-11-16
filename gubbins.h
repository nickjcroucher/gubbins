/*
 *  Wellcome Trust Sanger Institute
 *  Copyright (C) 2011  Wellcome Trust Sanger Institute
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

#ifndef _GUBBINS_H_
#define _GUBBINS_H_

#include "seqUtil.h"
#include "Newickform.h"

void run_gubbins(char vcf_filename[], char tree_filename[], char phylip_filename[]);
void extract_sequences(char vcf_filename[], char tree_filename[]);
char *calculate_ancestor_sequence(char * ancestor_sequence, char ** child_sequences, int sequence_length, int number_of_child_sequences);


#define MAX_EDGES_IN_TREE 10000

#endif



