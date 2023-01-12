// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

/**
 * @file tiff2ppm.c
 * 
 * \a tiff2ppm.exe - Convert a tiff image file into a ppm file (grey-scale, 16
 * bit file). Note that this conversion may lose information.
 * 
 * \par Usage: tiff2ppm.exe \<input tiff file\> \<output ppm file\> \par 
 *
 * \par Example:
 * \verbatim  tiff2ppm.exe my_reconstruction.tif my_reconstruction.ppm \endverbatim
 * 
 **/

#include <iostream>
#include <stdlib.h>
#include <Double_2D.h>
#include <io.h>

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //check for 3 arguements
  if(argc!=3 ){
    cout << "Wrong number of arguments. Usage: " 
	 << "tiff2ppm <input tiff file> <output ppm file>" << endl;
    return 1;
  }

  //read the data block in the file
  Double_2D data;
  int status;

  //read the data into an array
  status = read_tiff(argv[1], data);
  
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }
  
  //write the data to a file
  write_ppm(argv[2], data);
      
  return 0;
}
