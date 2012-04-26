#include <string>
#include <boost/python.hpp>

// don't include this
//#include "gubbins.h"
extern "C" void run_gubbins(char vcf_filename[], char tree_filename[], char phylip_filename[], char multi_fasta_filename[]);

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
		(char *) phylip_filename.c_str(),
		(char *) multi_fasta_filename.c_str());
  }
};

using namespace boost::python;
BOOST_PYTHON_MODULE(libPyGubbins) {
  class_<PyGubbins>("PyGubbins", init<>())
    .def("run", &PyGubbins::run)
    ;
}
