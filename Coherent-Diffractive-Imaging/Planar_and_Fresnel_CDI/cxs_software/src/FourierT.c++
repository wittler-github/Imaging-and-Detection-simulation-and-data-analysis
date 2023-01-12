#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "Complex_2D.h"
#include "FourierT.h"

using namespace std;


//constructor
FourierT::FourierT(int x_size, int y_size) 
  : nx(x_size),
    ny(y_size){

  //set up the fourier transforming array and plans
  original = new fftw_complex[nx*ny];
  transformed = new fftw_complex[nx*ny];

  //create the fftw plans. see FFTW3 documentation for an
  //explanation.
  f_forward = fftw_plan_dft_2d(nx, ny, original, transformed, 
			       FFTW_FORWARD, FFTW_MEASURE);
  f_backward = fftw_plan_dft_2d(nx, ny, transformed, original, 
				FFTW_BACKWARD, FFTW_MEASURE);
  
}

//clean up
FourierT::~FourierT(){

  fftw_destroy_plan(f_forward);
  fftw_destroy_plan(f_backward);

  delete[] original;
  delete[] transformed;
}


void FourierT::perform_forward_fft(Complex_2D & c_in){

  copy_to_fftw_array(original, c_in);
  fftw_execute(f_forward);
  copy_from_fftw_array(transformed, c_in);

}


void FourierT::perform_backward_fft(Complex_2D & c_in){

  copy_to_fftw_array(transformed, c_in); 
  fftw_execute(f_backward);
  copy_from_fftw_array(original, c_in);

}

void FourierT::copy_to_fftw_array(fftw_complex * array , Complex_2D & c){
  //check the dimensions:
  if(c.get_size_x()!=nx || c.get_size_y()!=ny ){
    cout << "In FourierT::copy_to_fftw_array. Dimensions of "
	 << "the Complex_2D and fftw_complex do not "
	 << "match.. exiting" <<endl;
    exit(1);
  }

  //quick copy
  memcpy(array,c.array,sizeof(fftw_complex)*nx*ny);

}

void FourierT::copy_from_fftw_array(fftw_complex * array, Complex_2D & c){
  //check the dimensions:
  //  if(c.get_size_x()!=nx || c.get_size_y()!=ny ){
  //  cout << "Dimensions of the Complex_2D and "
  //	 << "fftw_complex do not match.. exiting" <<endl;
    
  //    exit(1);
  //  }

  //always scale as FFTW doesn't normalise the result
  //we can't do a quick copy then :(

 double scale_factor = 1.0/(sqrt(nx*ny));
 
 for(int i=0; i < nx; ++i){
   for(int j=0; j < ny; ++j){
     c.set_real(i,j,(array[i*ny + j][REAL])*scale_factor);
     c.set_imag(i,j,(array[i*ny + j][IMAG])*scale_factor);
   }
 }
  
}

