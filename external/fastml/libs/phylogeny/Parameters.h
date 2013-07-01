#ifndef _Parameters_h
#define _Parameters_h

#include <iostream>
#include <ostream>
#include <string>
//#include "macros.h"
//#include "DebugStream.h"
//#include "StringUtils.h"

using std::string;
using std::istream;
using namespace std;

/*
CLASS
  Parameters

  A utility class used to manage program parameters. The class supports 
  setting default values for parameters, reading values from a parameters
  file and accessing parameters values from other parts of the program.

KEYWORDS
  parameters

AUTHORS
  Meir Fuchs (mailto: meirfux@math.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 1997.

CHANGES LOG
<UL>
<LI>9.01.05 Dina:
Bug fix: adding check to iterator end() to findInsertionPoint result
to paramType, getInt, getString, getFloat functions 
</LI>
<LI>17.05.04 Oranit Dror:
Adding new methods: dump() and empty()
</LI>
</UL>

GOALS
  Aid in managing program parameters. The Parameters class's main goal is to
  relieve programmers from the need to rewrite specialized parameters reading
  code sections for each of the programs. The Parameters class holds integer,
  floating point or string values in static storage indexed using the
  parameter's name. Class also supplies method for parsing strings.

USAGE
  The following section covers several issues regarding the Parameters class
  and its usage. Users should understand the issues covered below before
  using the class.

USAGE: SETTING DEFAULT PARAMETERS
  Default parameters are set using the addParameter methods. Note that the 
  type of the parameter is set according to the addParameter arguments. If
  a parameter is set using addParameter with an integer argument then
  subsequent updates (using updateParameter) to the same parameter will all 
  be stored as integers. Therefore the following code should output a 0:
  EXAMPLE
    Parameters::addParameter("Dummy", 3);
    Parameters::updateParameter("Dummy", "This should set it to zero");
    cout << Parameters::getstring("Dummy");
  END
  
  Note also that when setting defuault values of float parameters always use
  a decimal point or else these parameters will be added as intgers. For
  example:
  EXAMPLE
    Parameters::addParameter("CubeSize", 1.0);  OK
    Parameters::addParameter("CubeSize", 1);    Not OK. Integer parameter
  END

USAGE: READING PARAMETERS FROM FILE
  The readParameters method recieves an input stream from which parameters are
  to be read. Files are structured so that each line specifies the value of a
  parameter. Each line gives the parameter name, a white space and then the
  parameter value. Lines whose first non white-space charachter is # are
  ignored. A basic schema for using the Parameters class is to set the default
  values using addParameter calls and then calling readParameters to read in
  parameters with other values or new parameters. The following example works
  as such using the Parameters::dump method to print all the parameters
  and their values:
  EXAMPLE
    Parameters::addParameter("CubeSize", 1.0);
    Parameters::addParameter("MinVote", 8);
    ifstream params("params");
    Parameters::readParameters(params);
    params.close();
    Parameters::dump(cout);
  END
  With the following parameters file:
  EXAMPLE
    CubeSize 0.5
    File  pdb4hhb.ent
  END
  The following output should result:
  EXAMPLE
    CubeSize (Float) 0.5
    File     (Str)   pdb4hhb.ent
    MinVote  (Int)   8
  END

USAGE: ACCESSING PARAMETERS VALUES
  using the getInt, getFloat and getstring methods one may access the 
  parameters values. Note that a value will always be returned even if the 
  parameter is not stored as the same type. The get methods attempt to
  convert the parameter type to the requested return type of the method.
  The follwing code should produce 3 1's as its output:
  EXAMPLE:
    Parameters::addParameter("MaxMix", 1);   OK added an integer parameter
    cout << Parameters::getInt("MaxMix");
    cout << Parameters::getFloat("MaxMix");
    cout << Parameters::getstring("MaxMix");
  END
  Also note that parameters names are case sensitive.

USAGE: SUBCLASSING AND PERFORMANCE
  The Parameters engine keeps the parameters in a sorted list. Although
  finding a parameter and its value in this list is considerably fast most
  users will not want this overhead of searching for the parameter using
  string comparisons inside their main loops, as part of a code which can be 
  executed a great number of times. 
  The idea is to subclass the Parameters class and hold the values which
  require direct and fast access in seperate static variables. All parameters
  are accessed not throguh the getParameter methods but rather through
  specialized methods of the subclass. The following is an example of such an
  implementation. Notice the readParameters method.
  EXAMPLE:
    static int min_vote = 8;         // Default values
    static float cube_size = 1.0;

    class ProgParams : protected Parameters
    {
      int minVote()  { return min_vote };

      float cubeSize() { return cube_size };

      // file name is not held in static variable. Don't care about parameter
      // access time.
      string fileName() { return getstring("FileName"); }

      int readParameters(char* paramsfile) {
        addParameter("MinVote", min_vote);
        addParameter("CubeSize", cube_size);

        ifstream params(paramsfile);
        Parameters::readParameters(params);
        params.close();
        
        min_vote = getInt("MinVote");
        cube_size = getFloat("CubeSize");
      }
    } 
  END
*/
class Parameters 
{
public:
  //// Used by the paramType method. See below.
  enum ParamType { Undef, Int, Float, Str };

  //// readParameters recieves an input stream and reads parameters off this
  // input stream. See the usage section for details of how a parameters
  // file may be structured.
  static void readParameters(istream& paramStream);
  
  ////
  // Returns true if no parameters are defined. <br>
  // Author: Oranit Dror (oranit@tau.ac.il)   
  static bool empty(); 

  // GROUP: Setting parameters

  //// Adds an integer parameter. The integer value added will actually be
  // stored as an integer. Subsequent updates to the same parameter using
  // updateParameter will all be stored as integers.
  static void addParameter(const string& paramName, const int value);

  //// Adds a float parameter. The float value added will actually be
  // stored as a float. Subsequent updates to the same parameter using
  // updateParameter will all be stored as floats.
  static void addParameter(const string& paramName, const double value);

  //// Adds a string parameter. The string value added will actually be
  // stored as a string. Subsequent updates to the same parameter using
  // updateParameter will all be stored as strings.
  static void addParameter(const string& paramName, const string& value);

  //// Update the parameter value without changing the parameter type. The
  // value parameter is converted to the parameter's type if this parameter
  // already exists. If the parameter is not yet listed then updateParameter
  // adds a new parameter of string type.
  static void updateParameter(const string& paramName, 
                              const char* const value);

  // GROUP: Getting parameters values.
  
  //// Returns the storage type of the given parameter. If a parameter
  // of the given name does not exist then Undef is returned. See enum
  // ParamType above for possible return values.
  static ParamType paramType(const string& paramName);

  //// Gets the integer value of a given parameter. If parameter is not of
  // integer type then its value is converted to integer. If parameter does
  // not exist a 0 is returned.
  static int       getInt(const string& paramName, const int& defaultValue=0);

  //// Gets the float value of a given parameter. If parameter is not of
  // float type then its value is converted to float. If parameter does
  // not exist a 0 is returned.
  static float     getFloat(const string& paramName, const float& defaultValue=0.0);

  //// Gets the string value of a given parameter. If parameter is not of
  // string type then its value is converted to string. If parameter does
  // not exist an empty string is returned.
  static string    getString(const string& paramName, const string& defaultValue=string());

  // GROUP: Other methods



  //// Output all listed parameters. Used for debugging.
  static void   dump(ostream& out);

  //// Output all listed parameters. Used for debugging.
  //static void   dump(DebugStream& out, const unsigned int msgLevel);

  //// 
  // Output all listed parameters. <br>
  // Author: Oranit Dror (oranit@tau.ac.il) 
  static void dump(FILE* outputFile);

  //// A utility method. nextToken recieves an argument string, finds the first
  // white-space delimited token in this string and returns it while cutting 
  // this token off of the argument string (It it passed by reference). Tokens
  // are returned without any spaces. This method may be used repetitively to
  // tokenize a string.
  static string nextToken(string& str);

protected:
  //// Constructor is protected since all methods are static. No need to 
  // actually form an instance of this class.
  Parameters();
};

#endif


