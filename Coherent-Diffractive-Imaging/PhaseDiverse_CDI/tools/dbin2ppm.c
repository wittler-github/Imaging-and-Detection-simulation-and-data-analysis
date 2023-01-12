// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

//author: Nadia Davidson
//data: 17/01/2011

/**
 * @file dbin2ppm.c
 * 
 * \a dbin2ppm.exe - Convert a binary file (2D double / 64 bit, format) to
 * a ppm (grey-scale, 16 bit file). Note that this conversion loses
 * information.
 * 
 * \par Usage: dbin2ppm.exe \<input dbin file\> \<output ppm file\>
 * \<pixels in x\> \<pixels in y\> \par 
 *
 * \par Example:
 * \verbatim  dbin2ppm.exe my_reconstruction.dbin my_reconstruction.ppm 1024 1024 \endverbatim
 * 
 **/

#include <iostream>
#include <stdlib.h>
#include <Double_2D.h>
#include <io.h>

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //check for 5 arguements
  if(argc!=5 ){
    cout << "Wrong number of arguments. Usage: " 
	 << "dbin2ppm <input dbin file> <output ppm file> <dim x> <dim y>" << endl;
    return 1;
  }

  //read the data block in the file
  Double_2D data;
  int nx = atoi(argv[3]);
  int ny = atoi(argv[4]);
  int status;

  //read the data into an array

  status = read_dbin(argv[1], nx, ny, data);
  
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }
  
  //write the data to a file
  write_ppm(argv[2], data);
      
  return 0;
}
