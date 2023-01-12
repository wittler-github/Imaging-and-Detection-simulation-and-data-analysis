// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <cmath>
#include <ctype.h>
#include <io.h>
#include <Double_2D.h>

using namespace std;

#define FAILURE 0
#define SUCCESS 1

/***************************************************************/

/***************************************************************/
int write_spec(string file_name, const Double_2D & data){


  int nx = data.get_size_x();
  int ny = data.get_size_y();

  
   ofstream new_file;
   new_file.open(file_name.c_str());
   
   //error check.
   if(!new_file.is_open()){
     cout << "Could not open the file " 
       << file_name << "for writing" << endl;
     return FAILURE;
   }

   new_file << "Energy"<<'\t'<<"F.Density"<<endl; 

   for(int j=0; j < ny; ++j){
     for(int i=0; i < nx; ++i){

       new_file << data.get(i,j) <<'\t'; 

     }
     new_file << endl;
   }
   new_file.close();

   return SUCCESS; //success

}

/***************************************************************/

/***************************************************************/
int read_spec(string file_name, Double_2D & data){

  int ny=2;

  //open the input file:
  ifstream file(file_name.c_str());

  //error check.
  if(!file.is_open()){
    cout << "Could not open the file " << file_name << endl;
    return FAILURE;
  }

  string line, token;
  stringstream ssline;

  //make a temporary 1-d array;
  vector<string> string_data;

  while(getline(file,line)){

    ssline<<line;

    while(getline(ssline, token, '\t')){

      string_data.push_back(token);

    }

    ssline.clear();

    //only interested in lines which are numbers
    char c = (string_data.at(0)).at(0);

    if(string_data.size()>0 && !isdigit(c))
      string_data.clear();

  }

  file.close();

  int nx = string_data.size()/ny;

  if(nx>100){
    cout << "You are using " << nx << " number of points"
      <<"which is quite a lot. I usually find ~50 to be a "
      <<"good number."<< endl;
  }

  //do a sanity check

  //fill the output values
  if(data.get_size_x()==0)
    data.allocate_memory(nx,ny);

  float value = 0;

  for(int i=0,k=0; i< nx; ++i){
    for(int j=0; j < ny; ++j, ++k){
      std::istringstream(string_data.at(k)) >> value;
      //if(energy){
      if(j==0){
	value = 1.0/value;
      }
      // }
      data.set(i,j,value);
    }
  }

  return SUCCESS; //success

}

