#include <iostream>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "Config.h"

#define FAILURE 0
#define SUCCESS 1

#define KEY 0
#define VALUE 1

using namespace std;

Config::Config(string file_name){
  //open the file and
  //parse it's contents.
     
    ifstream file(file_name.c_str());
    if(!file.is_open()){
      cout << "Could not open the file " << file_name << endl;
      status = FAILURE;
      return;
    }

    string line;
    
    while(!file.eof()){
      getline(file,line);
      
      //tokenise the line of text

      //first remove the comments
      size_t pos = line.find("#");
      if(pos!=string::npos){
	line = line.substr(0,pos);
      }

      //now we get the key part and value part of the string
      pos = line.find("=");
      if(pos!=string::npos){
      string key_temp = line.substr(0,pos);

      //now get rid of white space      
      istringstream key_stream(key_temp);
      string key;      
      key_stream >> key;

      string value = line.substr(pos+1);

      mapping.insert(pair<string,string>(key,value));
      }
    }
    
    status=SUCCESS;
   
};
 


string Config::getString(string key){
  MapType::iterator iter = mapping.find(key);
  if (iter == mapping.end() ){
    std::cout << "Key:"<<key<<" not found in config file" << endl;
    status=FAILURE;
    return "";
  }

  //remove white space
  string value;      
  istringstream value_temp(iter->second);
  value_temp >> value;
  return value;
}

list<string> * Config::getStringList(string key){
  MapType::iterator iter = mapping.find(key);
  if (iter == mapping.end() ){
    std::cout << "Key:"<<key<<" not found in config file" << endl;
    status=FAILURE;
    return 0;
  }
  //tokenize the string
  list<string> * new_list = new list<string>;
  string value;
  istringstream value_temp(iter->second);
  value_temp >> value;
  while(value_temp){
    new_list->push_back(value);
    value_temp >> value;
  }
  if(new_list->size()==0){
    std::cout << "Key:"<<key
	      <<" has no entries in the config file" << endl;
    status=FAILURE;
    return 0;
  }
  
  return new_list;
}

list<int> * Config::getIntList(string key){
  MapType::iterator iter = mapping.find(key);
  if (iter == mapping.end() ){
    std::cout << "Key:"<<key<<" not found in config file" << endl;
    status=FAILURE;
    return 0;
  }
  //tokenize the string
  list<int> * new_list = new list<int>;
  string value;
  istringstream value_temp(iter->second);
  value_temp >> value;
  while(value_temp){
    new_list->push_back(atoi(value.c_str()));
    value_temp >> value;
  }
  if(new_list->size()==0){
    std::cout << "Key:"<<key<<" has no entries in the config file" << endl;
    status=FAILURE;
    return 0;
  }
  return new_list;
}

double Config::getDouble(string key){
MapType::iterator iter = mapping.find(key);
 if (iter == mapping.end() ){
   std::cout << "Key:"<<key<<" not found in config file" << endl;
   status=FAILURE;
   return 0;
 }
 //convert to a double
 return atof((iter->second).c_str());
 
}

int Config::getInt(string key){
MapType::iterator iter = mapping.find(key);
 if (iter == mapping.end() ){
   std::cout << "Key:"<<key<<" not found in config file" << endl;
   status=FAILURE;
   return 0;
 }
 //convert to an int
 return atoi((iter->second).c_str());
}

