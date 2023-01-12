/**
 * @file FresnelCDI_simulation_example.c
 *
 * This file is a simulation of fCDI for the purpose of imaging small 
 * contrast changes in biological samples. It is based on Nadia's example 
 * and uses her header files.
 * see: FresnelCDI_example.c
 *
 * @author Michael Jones <michael.jones@latrobe.edu.au>,
 * Nadia Davidson <nadiamd@unimelb.edu.au> 
 *
 * Last modified on 18/3/2011
 *
 */

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "io.h"
//#include "utils.h"
#include "Complex_2D.h"
#include "Double_2D.h"
#include "FresnelCDI.h"
#include <cstdlib>
#include <math.h>
#include <string>

using namespace std;

/**********************************/

int main(int argc, char * argv[]){

  /** All of this part is just reading files and setting constants **/

  //number of output iterations
  int output_iterations = 25;
  int total_iterations = 200;

  //load an image file for the sample
  Double_2D object;
  read_image("image_files/FCDI_simulation_object.tiff", object);

  //get the object dimensions
  int nx = object.get_size_x();
  int ny = object.get_size_y();

  //load the support to use in the reconstruction
  Double_2D support;
  read_image("image_files/FCDI_simulation_support.tiff", support);
  
  //we will use the same white field as the Fresnel examples. For
  //this reason you will need to reconstruct it first using
  //FresnelCDI_WF_example.exe
  Complex_2D wf(nx,ny);
  int status = read_cplx("wf_recovered.cplx", wf);
  if(!status){ //give an error if we couldn't open the file.
    cout  << "Maybe you need to run ./FCDI_WF_example.exe "
	  << "first... exiting"  << endl;
    return(1);
  }
  
  //make sure the object and support are the same dimensions
  if(support.get_size_x() != nx || support.get_size_y() != ny){
    cout << "dimensions do not match....exiting" <<endl;
    return(1);
  }


  /**********Set up the projection*************************
   This is like creating a virtual experiment. Set the
   experimental parameters etc. and get it ready to simulate
   the diffration pattern (and later we will use it for
   reconstruction as well).
  *********************************************************/

  //make a 2D array and allocate some memory.
  //This will be used to output the image of the 
  //current exit-surface-wave estimate.
  Double_2D result(nx,ny);
  
  //set the experimental parameters (all are in meters)
  double wavelength = 4.892e-10; //wavelength
  double fd = 0.8932; //focal to detector
  double fs = 2.16e-3; //focal to sample
  double ps = 13.5e-6; //pixel size

  //create the projection object which will be used later to
  //perform the reconstuction.
  Complex_2D object_estimate(nx,ny);

  FresnelCDI proj(object_estimate, //estimate
		  wf, //white field 
		  wavelength, //wavelength
		  fd, //focal-to-detector
		  fs, //focal-to-sample
		  ps, //pixel size
		  1.0); //normalisation of white-field to sample data

  //***************************************************
  //  Simulate a diffraction pattern
  //***************************************************  

  /*** First we will create a transmission function **/

  // create the complex array which it will be stored in 
  Complex_2D input(ny,ny);

  //lets say our sample is 150nm maximum thickness and is made from gold
  double thickness = 150e-9;
  double delta = 6.45e-4;
  double beta = 1.43e-4;
  double k = (2.0 * M_PI)/wavelength;

  double max = object.get_max();

  for (int i = 0; i<nx; i++){
    for (int j = 0; j<ny; j++){

      //calculate the amplitude and phase of the transmission function
      //at each point (from the image).
      double this_pixel_thickness = thickness*object.get(i,j)/max;
      double mag = exp(-beta*k*this_pixel_thickness);
      double phase = -delta*k*this_pixel_thickness;

      input.set_value(i,j,MAG, mag);
      input.set_value(i,j,PHASE, phase);
    }
  }

  //save the amplitude and phase to file so we know how it looks
  input.get_2d(MAG,result);
  write_tiff("trans_mag_input.tiff",result);

  input.get_2d(PHASE,result);
  write_tiff("trans_phase_input.tiff",result);

  /********* now simulate: ***********************************/ 
   
  //get the illuminating white field in the sample plane
  proj.propagate_from_detector(wf);
  
  //multiply the transmission function by the white field
  input.multiply(wf);

  //get the magnitude of the wave
  input.get_2d(MAG_SQ,result);
  //write the output to file
  write_tiff("sample_plane.tiff",result);

  //propagate to the detector
  proj.propagate_to_detector(input);

  //take the wf back to the detector so we're ready for the reconstructions
  proj.propagate_to_detector(wf);

  //write input intensity to Double_2D array
  //"result" will be used as the diffraction data.
  input.get_2d(MAG_SQ,result);

  //write the diffraction pattern to a file
  write_tiff("forward_projection.tiff",result);

  //if you feel like making a more realistic simulation with threasholding
  //and simulated noise, consider uncommenting the code below.
  //add random noise to the output. Where, and in what form is best??
  /**for(int i=0; i<nx; i++){
     for(int j=0; j<ny; j++){
        result.set(i,j,result.get(i,j) * (rand()%21) );
	if(result.get(i,j)<0)
      	   result.set(i,j,0);
     }
     }**/
  
  //the approx. theshold level of the image.
  //Not too much (less than 1!, depending on absorption level below)
  //const double noise_level = 0.5;
  
  //apply a threshold to make the simulation a bit more realistic
  /**  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      result.set(i,j,result.get(i,j)-noise_level);
      if(result.get(i,j)<0)
	result.set(i,j,0);
	}
	}**/
  
  
  /****************************************************/
  /*******  set up the reconstuction ****************
   **using the diffraction pattern produced by the previous 
   **projection as the starting point***********************/
  /******************************************************/
  /*******************************************************/


  //set the support and intensity
  proj.set_support(support);  
  proj.set_intensity(result);
  
  //set the algorithm
  proj.set_algorithm(ER);
  
  //Initialise the current object ESW
  proj.initialise_estimate(0);
  
  //set a complex contraint for the transmission function
  //proj.set_complex_contraint_function(charge_flip);
  
  /*** run the reconstruction ************/
  
  for(int i=0; i<  total_iterations+1; i++){

    cout << "iteration " << i << endl;

    //apply the iterations  
    proj.iterate(); 
    cout << "Error: " << proj.get_error() << endl;

    if(i%output_iterations==0){

      //output the current estimate of the object
      object_estimate.get_2d(MAG,result);
      ostringstream temp_str ( ostringstream::out ) ;
      temp_str << "fcdi_example_iter_" << i << ".tiff";
      write_tiff(temp_str.str(),result);
    
      //uncomment to apply the shrinkwrap algorithm
      //      proj.apply_shrinkwrap();
    }
  }

  //get the reconstructed transmission function for 
  //the final iteration
  Complex_2D trans(nx,ny);
  proj.get_transmission_function(trans);

  //and write the magnitude and phase of it to an image file
  trans.get_2d(MAG,result);
  write_tiff("trans_mag_recovered.tiff",result);

  trans.get_2d(PHASE,result);
  write_tiff("trans_phase_recovered.tiff",result);

  return 0;
  
}
