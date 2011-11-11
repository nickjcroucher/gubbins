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

#ifndef _SNP_SITES_H_
#define _SNP_SITES_H_

int build_reference_sequence(char reference_sequence[], FILE * alignment_file_pointer);
int detect_snps(char reference_sequence[], FILE * alignment_file_pointer, int length_of_genome);
void build_snp_locations(int snp_locations[], char reference_sequence[]);
void get_bases_for_each_snp(FILE * alignment_file_pointer, int snp_locations[], char ** bases_for_snps, int length_of_genome, int number_of_snps);
int generate_snp_sites(char filename[]);

#endif