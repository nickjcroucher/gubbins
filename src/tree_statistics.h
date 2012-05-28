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

#ifndef _TREE_STATISTICS_H_
#define _TREE_STATISTICS_H_

 void create_tree_statistics_file(char filename[], sample_statistics ** statistics_for_samples, int number_of_samples);
 float recombination_to_mutation_ratio(int number_of_recombinations, int number_of_snps);
 float recombination_blocks_to_mutation_ratio(int number_of_blocks, int number_of_snps);
#define MAX_FILE_NAME_SIZE 1024

#endif




