#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <cmath>
#include "io.h"
#include "Complex_2D.h"
#include <math.h>

using namespace std;

#define FAILURE 0
#define SUCCESS 1

/***************************************************************/

/***************************************************************/
int read_cplx(string file_name, Complex_2D & complex){

  int nx = complex.get_size_x();
  int ny = complex.get_size_y();

  //open the input file:
  FILE * file = fopen(file_name.c_str(), "r+b");

  //error check.
  if(!file){
    cout << "Could not open the file " << file_name << endl;
    return FAILURE;
  }

  double * buffer = new double [nx*ny*2];
  size_t elements_read = fread(buffer, sizeof(double), nx*ny*2, file);

  if(elements_read!=(nx*ny*2)){
    cout << "Could not correctly read the file " << file_name
	 << ". Perhaps the dimensions specifies are wrong?"<<endl;
    return FAILURE;
  }

  fclose(file);
  
  //Do a sanity check. Is the file size right 
  //for a nx by ny array with 16 bit pixel values?
  
  for(int i=0; i < nx; ++i){
    for(int j=0; j< ny; ++j){
      complex.set_real(i,j,buffer[2*((ny-j-1)*nx+i)]);
      complex.set_imag(i,j,buffer[2*((ny-j-1)*nx+i)+1]);
    }
  }
  
  delete buffer;

  return SUCCESS; //success
    
}

/***************************************************************/

/***************************************************************/
int write_cplx(string file_name, const Complex_2D & complex){

  int nx = complex.get_size_x();
  int ny = complex.get_size_y();

  //open the input file:
  FILE * file = fopen(file_name.c_str(), "w+b");

  //error check.
  if(!file){
    cout << "Could not open the file " << file_name << endl;
    return FAILURE;
  }

  double * buffer = new double [nx*ny*2];

  //Do a sanity check. Is the file size right 
  //for a nx by ny array with 16 bit pixel values?
  
  for(int i=0; i < nx; ++i){
    for(int j=0; j< ny; ++j){
      buffer[2*((ny-j-1)*nx+i)] = complex.get_real(i,j);
      buffer[2*((ny-j-1)*nx+i)+1] = complex.get_imag(i,j);
    }
  }
  
  fwrite(buffer, sizeof(double), nx*ny*2, file);
  fclose(file);
  
  delete buffer;

  return SUCCESS; //success
}
