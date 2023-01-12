// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

/**
 * @file Config.h
 * @class Config
 * @author Nadia Davidson
 * @date Last Modified on 6/1/2001
 *
 * @brief A configuration file parser
 *
 * This class allows simple text files to be read and is used by the
 * command line tool programs.  The parser takes a file name and
 * constructs key-value pairs for each line with the syntax: 
 *
 * key = value   or
 *
 * key = value_1 value_2 value_3 ...
 *
 * White spaces are ignored and everything to the right of a "#"
 * symbol is also ignored (this allows comments to be included in the
 * config file. See example.config to see how it looks.  If a
 * key (as given by any of the methods below) is not found in the
 * file, the status is set to FAILURE and the method returns either 0,
 * an empty string or list.
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <map>
#include <list>

#define FAILURE 0
#define SUCCESS 1

typedef std::map<std::string,std::string> MapType;

class Config{

 private:

  /** The key-value map */
  MapType mapping;

  /** Were all the requested fields found in the file? */
  int status;

 public:

  /**
   * Construct a configuration object.
   * 
   * @param file_name The name of the file to parse.
   */
  Config(std::string file_name);

  /**
   * Given a key, return the value from the config file as a string
   * type.
   * 
   * @param key The key, or field name
   * @return The value 
   */
  std::string getString(std::string key);

  /**
   * Given a key, return the value from the config file as a double
   * type.
   * 
   * @param key The key, or field name 
   * @return The value
   */
  double getDouble(std::string key);

  /**
   * Given a key, return the value from the config file as an integer
   * type.
   * 
   * @param key The key, or field name 
   * @return The value
   */
  int getInt(std::string key);

  /**
   * Given a key, return a list of values from the config file as
   * integer types.
   * 
   * @param key The key, or field name 
   * @return The value
   */
  std::list<int> * getIntList(std::string key);

  /**
   * Given a key, return a list of values from the config file as
   * string types.
   * 
   * @param key The key, or field name 
   * @return The value
   */  
  std::list<std::string> * getStringList(std::string key);

  /**
   * Were all the requests for a key-value pair successful?  If any
   * of the keys were not found in the file (after calling the access
   * methods above), the status is set to FAILURE (0).
   * 
   * @return Whether all requests were successful.
   */
  int getStatus(){return status;};

    
};

#endif
