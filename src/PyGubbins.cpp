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
#include <string>
#include <boost/python.hpp>

// don't include this
//#include "gubbins.h"
extern "C" void run_gubbins(char vcf_filename[], char tree_filename[], char multi_fasta_filename[], int min_snps);

class PyGubbins {
public:
  PyGubbins () {}; // empty default constructor
  void run(const std::string & vcf_filename,
	   const std::string & tree_filename, 
	   const std::string & phylip_filename, 
	   const std::string & multi_fasta_filename) {
    // todo: make run_gubbins take const char *
    run_gubbins((char *) vcf_filename.c_str(),
		(char *) tree_filename.c_str(),
		(char *) multi_fasta_filename.c_str()
		);
  }
};

using namespace boost::python;
BOOST_PYTHON_MODULE(libPyGubbins) {
  class_<PyGubbins>("PyGubbins", init<>())
    .def("run", &PyGubbins::run)
    ;
}
