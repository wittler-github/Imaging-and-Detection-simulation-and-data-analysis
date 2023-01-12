/**
 * @file hdf2ppm.c
 * 
 * \a hdf2ppm.exe - Extract the data from a HDF4 file and save as 
 * a ppm (grey-scale, 16 bit file). 
 * 
 * \par Usage: hdf2ppm.exe \<input hdf file\> \<output ppm file\> \<optional: data block name\> \par 
 *
 * \par Example:
 * \verbatim  hdf2ppm.exe my_data.hdf my_data.ppm  \endverbatim
 * 
 **/

#include <iostream>
#include <stdlib.h>
#include "Double_2D.h"
#include "io.h"

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //check for 2 or 3 arguements
  if(argc<3 || argc>4 ){
    cout << "Wrong number of arguments. Usage: " 
	 << "hdf2ppm <input hdf file> <output ppm file>" << endl;
    return 1;
  }

  //read the data block in the file
  Double_2D data;
  int status;

  //read the data into an array
  if(argc==4)
    status = read_hdf4(argv[1], data, argv[3]);
  else
    status = read_hdf4(argv[1], data);
  
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }
  
  //write the data to a file
  write_ppm(argv[2], data);
      
  return 0;
}
