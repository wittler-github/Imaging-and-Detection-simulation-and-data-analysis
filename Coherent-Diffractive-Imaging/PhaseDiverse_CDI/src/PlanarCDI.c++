// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

#include <cstdlib> 
#include <cmath>
#include <Complex_2D.h>
#include <PlanarCDI.h>

//double ** PlanarCDI::get_intensity_autocorrelation(){
void PlanarCDI::get_intensity_autocorrelation(Double_2D & autoc){

  //make a temporary object
  Complex_2D temp_intensity(nx,ny);

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      //set the real and imaginary components using the magnitude.
      double component = (1.0/sqrt(2.0))*(intensity_sqrt.get(i,j)*intensity_sqrt.get(i,j));
      temp_intensity.set_value(i,j,REAL, component);
      temp_intensity.set_value(i,j,IMAG, component);
    }
  }
  
  // fourier transform the intensity 
  //propagate_from_detector(temp_intensity);  

  temp_intensity.perform_backward_fft(); 
  temp_intensity.invert(true);

  //get the magnitude of the fourier transformed data.
  temp_intensity.get_2d(MAG, autoc);
  
}


void PlanarCDI::initialise_estimate(int seed){
  //initialise the random number generator
  srand(seed);

  int max_value = intensity_sqrt.get_max();

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      if(!support.get(i,j)){ //enforce the support condition on the inital guess
	complex.set_value(i,j,REAL,0); 
	complex.set_value(i,j,IMAG,0);
      }
      else{
	double r = support.get(i,j)*(max_value*rand()/(double) RAND_MAX);
	double im = support.get(i,j)*(max_value*rand()/(double) RAND_MAX);

	complex.set_value(i,j,REAL,r); 
	complex.set_value(i,j,IMAG,im);
      }
    }
  }

}


void PlanarCDI::propagate_to_detector(Complex_2D & c){
    c.perform_forward_fft();
    c.invert(true);
}

void PlanarCDI::propagate_from_detector(Complex_2D & c){
  c.invert(true);
  c.perform_backward_fft(); 
}
