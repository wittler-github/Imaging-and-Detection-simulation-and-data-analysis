/**
 * @file cplx2dbin.c
 * 
 * \a cplx2dbin.exe - Extract part of a complex binary file (2D fftw format)
 * and save as a real binary file (2D double / 64 bit, format). The real, 
 * imaginary, magnitude, phase and magnitude squared can be extracted.
 * 
 * \par Usage: cplx2dbin.exe \<input cplx file\> \<output dbin file\> \<component type\> \<size in x\> \<size in y\> 
 * \par
 * where component type is one of:
 * - 0: REAL 
 * - 1: IMAG 
 * - 2: MAG 
 * - 3: PHASE 
 * - 4: MAG_SQ 
 *
 * \par Example:
 * \verbatim cplx2dbin.exe my_white_field.cplx my_white_field_illum.bin 4 1024 1024 \endverbatim
 * Extract the magnitude squared from a reconstucted white-field file. 
 * 
 **/





#include <iostream>
#include <stdlib.h>
#include "io.h"
#include "Complex_2D.h"
#include "Double_2D.h"

using namespace std;

/**************************************/
int main(int argc, char * argv[]){

  //check for 5 arguements
  if(argc!=6 ){
    cout << "Wrong number of arguments." <<endl;
    cout << "Usage: cplx2dbin <input cplx file> <output dbin file> "
	 << "<component type> <dim x> <dim y>" << endl;
    cout << "  where component type is one of:" <<endl;
    cout << "     0 - REAL" << endl;
    cout << "     1 - IMAG" << endl;
    cout << "     2 - MAG" << endl;
    cout << "     3 - PHASE" << endl;
    cout << "     4 - MAG_SQ" << endl;
    return 1;
  }

  //read the array size from the command arguments
  int nx = atoi(argv[4]);
  int ny = atoi(argv[5]);

  //do some initalisation
  Complex_2D data(nx,ny);
  int status;

  //read the data block into an array
  status = read_cplx(argv[1], data);
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }
  
  //Extract the required part of the complex field
  Double_2D real(nx,ny);
  data.get_2d(atoi(argv[3]),real);

  //write the result to a new file
  write_dbin(argv[2], real);
      
  return 0;
}
