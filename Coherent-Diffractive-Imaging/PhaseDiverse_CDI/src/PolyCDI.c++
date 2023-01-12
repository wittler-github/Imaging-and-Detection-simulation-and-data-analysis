// Copyright 2012 T'Mir Danger Julius for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include <stdlib.h>
#include <cstdlib> 
#include <Complex_2D.h>
#include <Double_2D.h>
#include <PolyCDI.h>
#include <io.h> 
#include <sstream>
#include <typeinfo>
#include <utils.h>

using namespace std;

//constructor for the class which handles partially coherent diffractive imaging.

PolyCDI::PolyCDI(Complex_2D & initial_guess,
    double beta, 
    int n_best,
    bool parallel
    )
:BaseCDI(initial_guess,n_best),
  beta(beta), 
  parallel(parallel),
  intensity_sqrt_calc(nx, ny){

  }

//destructor for cleaning up
PolyCDI::~PolyCDI(){
}

//Initialise the spectrum Double_2D from a file
void PolyCDI::set_spectrum(string file_name){

  read_spec(file_name, spectrum);
  nlambda=spectrum.get_size_x();

  double maxweight=0;
  double totweight=0;

  for(int i=0; i<nlambda; i++){
    if(spectrum.get(i, WEIGHT) > maxweight){
      lambdac = spectrum.get(i, WL);
      maxweight=spectrum.get(i, WEIGHT);
    }
    totweight+=spectrum.get(i, WEIGHT);
  }

  //std::cout<<"lambdac is "<<lambdac<<endl;

  for(int i=0; i<nlambda; i++){
    spectrum.set(i, WEIGHT, spectrum.get(i, WEIGHT)/totweight);
  }

}

//Initialise the spectrum Double_2D from a Double_2D
void PolyCDI::set_spectrum(Double_2D spectrum_in){
  spectrum=spectrum_in;
  nlambda=spectrum.get_size_x();

  double maxweight=0;
  double totweight=0;

  for(int i=0; i<nlambda; i++){
    if(spectrum.get(i, WEIGHT) > maxweight){
      lambdac = spectrum.get(i, WL);
      maxweight=spectrum.get(i, WEIGHT);
    }
    totweight+=spectrum.get(i, WEIGHT);
  }

  //std::cout<<"lambdac is "<<lambdac<<endl;

  for(int i=0; i<nlambda; i++){
    spectrum.set(i, WEIGHT, spectrum.get(i, WEIGHT)/totweight);
  }

}


//set the initial guess.
//fills the transmission function with a random 
//value between 0 and 1 within the support, and
//0 everywhere else. Complex is then set to be 
//the value of the transmission function
void PolyCDI::initialise_estimate(int seed){

  //initialise the random number generator
  srand(seed);

  int max_value = intensity_sqrt.get_max();

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      if(!support.get(i,j)){ 

	complex.set_value(i,j,REAL,0);
	complex.set_value(i,j,IMAG,0);
      }
      else{
	double re = support.get(i,j)*(max_value*rand()/(double) RAND_MAX);
	double im = support.get(i,j)*(max_value*rand()/(double) RAND_MAX);

	complex.set_value(i,j,REAL,re);
	complex.set_value(i,j,IMAG,im);
      }
    }
  }
  //  complex=transmission;
}

//this iterate function overrides that of BaseCDI, 
//and handles the multiple modes.
int PolyCDI::iterate(){

  //below is the code for the special case of ER
  //this is faster than using the generic algorithm code
  //further down in this function.

  Complex_2D temp_complex_PFS(nx, ny);
  Complex_2D temp_complex_PF(nx, ny);
  Complex_2D temp_complex_PS(nx, ny);
  Complex_2D temp_complex_PSF(nx, ny);

  singleCDI.clear();


  if(algorithm==ER){

    propagate_to_detector(complex);

    scale_intensity(complex);

    propagate_from_detector(complex);
    apply_support(complex);
    update_n_best();

    return SUCCESS;

  }

  //start of the generic algorithm code

  //PFS
  if(algorithm_structure[PFS]!=0){
    temp_complex_PFS=complex;
    apply_support(temp_complex_PFS);
    propagate_to_detector(temp_complex_PFS);
    scale_intensity(temp_complex_PFS);
    propagate_from_detector(temp_complex_PFS);
  }

  //F
  if(algorithm_structure[PF]!=0){
    temp_complex_PF=complex;
    propagate_to_detector(temp_complex_PF);
    scale_intensity(temp_complex_PF);
    propagate_from_detector(temp_complex_PF);
  }

  //S
  if(algorithm_structure[PS]!=0){
    temp_complex_PS=complex;
    apply_support(temp_complex_PS);
  }

  //SF
  if(algorithm_structure[PSF]!=0){
    if(algorithm_structure[PF]!=0){
      temp_complex_PSF=temp_complex_PF;
      apply_support(temp_complex_PSF);
    }else{
      temp_complex_PSF=complex;
      propagate_to_detector(temp_complex_PSF);
      scale_intensity(temp_complex_PSF);
      propagate_from_detector(temp_complex_PSF);
      apply_support(temp_complex_PSF);
    }
  }

  //combine the result of the separate operators
  double value_real, value_imag;
  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){

      //Add the identity
      value_real = (1+algorithm_structure[PI])*complex.get_real(i,j);
      value_imag = (1+algorithm_structure[PI])*complex.get_imag(i,j);

      //Add the component from the PfPs operator
      if(algorithm_structure[PFS]!=0){
	value_real+=algorithm_structure[PFS]*temp_complex_PFS.get_real(i,j);
	value_imag+=algorithm_structure[PFS]*temp_complex_PFS.get_imag(i,j);
      }

      //Add the component from the Pf operator
      if(algorithm_structure[PF]!=0){
	value_real+=algorithm_structure[PF]*temp_complex_PF.get_real(i,j);
	value_imag+=algorithm_structure[PF]*temp_complex_PF.get_imag(i,j);

      }

      //Add the component from the Ps operator
      if(algorithm_structure[PS]!=0){
	value_real+=algorithm_structure[PS]*temp_complex_PS.get_real(i,j);
	value_imag+=algorithm_structure[PS]*temp_complex_PS.get_imag(i,j);
      }

      //Add the component from the PsPf operator
      if(algorithm_structure[PSF]!=0){
	value_real+=algorithm_structure[PSF]*temp_complex_PSF.get_real(i,j);
	value_imag+=algorithm_structure[PSF]*temp_complex_PSF.get_imag(i,j);
      }

      complex.set_real(i,j,value_real);
      complex.set_imag(i,j,value_imag);

    }
  }

  //Update the transmission using the dominant mode
  update_n_best();
  return SUCCESS;
}

//scale the highest occupancy mode 
//this overwrites the function of the same
//name in BaseCDI
void PolyCDI::scale_intensity(Complex_2D & c){
  double norm2_mag=0;
  double norm2_diff=0;
  double measured_int=0;
  double calc_int=0;
  double scaled_mag=0;

  //  double lambdmax = spectrum.get(nlambda, wl);
  //  double nxmin = (1.0-lambdac/lambdmax)*((nx-1)/2);
  //  double nymin = (1.0-lambdac/lambdmax)*((ny-1)/2);

  expand_wl(c);

  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){
      //reset the magnitude
      if(beam_stop==0 || beam_stop->get(i,j)>0){

	calc_int=intensity_sqrt_calc.get(i,j);
	measured_int=intensity_sqrt.get(i,j); 

	if(calc_int==0.0){
	  scaled_mag=0.0;
	}else{
	  scaled_mag=c.get_mag(i,j)*measured_int/calc_int;
	}
	c.set_mag(i,j,scaled_mag);

	//calculate the error
	norm2_mag += measured_int*measured_int;

	norm2_diff += (measured_int-calc_int)*(measured_int-calc_int);

      }
    }
  }
  current_error = norm2_diff/norm2_mag;
}

//take the central wavelength, and generate the
//diffraction pattern 
void PolyCDI::expand_wl(Complex_2D & c){

  Double_2D intensity(nx, ny);

  /*  for(int i=0; i<nx; i++){
      for(int j=0; j<ny; j++){

      intensity.set(i, j, 0.0);
      }
      }
   */
  for(int sn=0; sn<nlambda; sn++){

    double lambda=spectrum.get(sn, WL);
    double weightl=spectrum.get(sn, WEIGHT);

    double lf=/*(2*lambdac-lambda)*/lambda/lambdac;

    for(int i=0; i<nx; i++){
      for(int j=0; j<ny; j++){

	double xval = lf*(i-nx/2.0)+nx/2.0;
	double yval = lf*(j-ny/2.0)+ny/2.0;

	if(((xval>=0) && (xval<(nx)))&&((yval>=0)&&(yval<(ny)))){

	  intensity.set(i, j, intensity.get(i, j) + weightl*pow(c.get_mag(xval, yval), 2));

	}
      }
    }
  }
  for(int i=0; i < nx; i++){
    for(int j=0; j < ny; j++){
      intensity.set(i, j, sqrt(intensity.get(i, j)));
    }
  }

  intensity_sqrt_calc=intensity;
}




//add the intensities across all modes scaled by the eigenvalues
Double_2D PolyCDI::sum_intensity(vector<Complex_2D> & c){

  Double_2D magnitude(nx, ny);

  //set the intensities to 0
  /*  for(int i=0; i< nx; i++){
      for(int j=0; j< ny; j++){
      magnitude.set(i,j,0.0);
      }
      }*/


  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){

      double mag_total=0;
      for(int mode=0; mode<c.size(); mode++){

	mag_total += eigen.at(mode)*c.at(mode).get_mag(i,j)*c.at(mode).get_mag(i,j);

      }
      magnitude.set(i,j,mag_total);
    }
  }

  return(magnitude);
};


//Propagate from the object plane to the detector
void PolyCDI::propagate_to_detector(Complex_2D & c){

  c.perform_forward_fft();
  c.invert(true);

}

//Propagate form the detector plane to the object
void PolyCDI::propagate_from_detector(Complex_2D & c){

  c.invert(true);
  c.perform_backward_fft();

}
