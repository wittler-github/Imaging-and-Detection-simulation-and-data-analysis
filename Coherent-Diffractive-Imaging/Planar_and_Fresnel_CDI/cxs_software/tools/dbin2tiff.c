//author: Nadia Davidson
//data: 23/02/2011

/**
 * @file dbin2tiff.c
 * 
 * \a dbin2tiff.exe - Convert a binary file (2D double / 64 bit, format) to
 * a tiff (grey-scale, 16 bit file).
 * 
 * \par Usage: dbin2tiff.exe \<input dbin file\> \<output tiff file\>
 * \<pixels in x\> \<pixels in y\> \par 
 *
 * \par Example:
 * \verbatim  dbin2tiff.exe my_reconstruction.dbin my_reconstruction.tiff 1024 1024 \endverbatim
 * 
 **/

#include <iostream>
#include <stdlib.h>
#include "Double_2D.h"
#include "io.h"

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //check for 5 arguements
  if(argc!=5 ){
    cout << "Wrong number of arguments. Usage: " 
	 << "dbin2tiff <input dbin file> <output tiff file> <dim x> <dim y>" << endl;
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
  status = write_tiff(argv[2], data);
      
  return 0;
}
