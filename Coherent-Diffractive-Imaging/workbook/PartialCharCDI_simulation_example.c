// Copyright 2011 Nadia Davidson 
// for The ARC Centre of Excellence in Coherent X-ray Science. 
//
// This program is distributed under the GNU General Public License. 
// We also ask that you cite this software in publications where you made 
// use of it for any part of the data analysis.
//
// date last modified: 17/01/2011

/**
 * @file PlanarCDI_simulation_example.c
 *
 * \a PlanarCDI_simulation_example.c The object from the
 * PlanarCDI_example (in fact I used the reconstructed image as the
 * object) is used to simulate a diffraction pattern for planar
 * CDI. The diffraction pattern is thresholded to make it more
 * realistic and then CDI reconstruction is performed.
 *
 */

#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib> 
#include <io.h>
#include <utils.h>
#include <Complex_2D.h>
#include <Double_2D.h>
#include <PartialCharCDI.h>
#include <sstream>

using namespace std;


/**************************************/
int main(void){

  //define some constants which will be used in the code:

  //the data file name
  const static char * data_file_name = "image_files/object.tiff";

  //the approx. level of background noise in the image
  const double noise_level = 15;

  //the file which provides the support (pixels with the value 0
  //are considered as outside the object)
  const static char * support_file_name = "image_files/planar_support.tiff";

  //number of hybrid input-out iterations to perform.
  const int hio_iterations = 0;

  //the number of error-reduction iterations to perform.
  const int er_iterations = 0;

  //output the current image ever "output_iterations"
  const int output_iterations = 10;

  //Coherence lengths of simulated beam (in meters)
  double lx = 2.0e-5;
  double ly = 2.0e-5;
  
  //Pixel size detector in m
  double psize_x=13.5e-6;
  double psize_y=13.5e-6;

  //Energy of the beam in eV
  double e_beam=1400;

  //Distance between detector and sample in metres
  double z_sd=1.4;

  
  /****** get the object from an image file ****************/

  //get the data from file
  Double_2D data;

  //read the data into an array
  int status = read_tiff(data_file_name, data);
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }

  int n_x = data.get_size_x();
  int n_y = data.get_size_y();
  
  /**** create the projection/reconstruction object *****/

  Complex_2D first_guess(n_x,n_y);
  PartialCharCDI my_partial(first_guess, z_sd, e_beam, psize_x, psize_y);

  // Convert coherence lengths in object-plane-meters to detector-pixels (using the native conversion methods in PartialCharCDI):
  my_partial.set_initial_coherence_guess_in_m(lx, ly);
  lx = my_partial.get_x_coherence_length_in_pixels();
  ly = my_partial.get_y_coherence_length_in_pixels();

  data = gaussian_convolution(data, lx, ly); // Apply simulated incoherence-blur to the object before feeding it to the CDI algorithm

  //fill the complex no. with image data
  Complex_2D input(n_x,n_y);
  for(int i=0; i<n_x; i++){
    for(int j=0; j<n_y; j++){
      input.set_value(i,j,REAL, 1.0/sqrt(2.0)*data.get(i,j));
      input.set_value(i,j,IMAG, 1.0/sqrt(2.0)*data.get(i,j));
    }
  }

  //propagate to the detector plane
  my_partial.propagate_to_detector(input);
  

  //write the fourier transform to file.
  Double_2D intensity(n_x,n_y);
  input.get_2d(MAG_SQ,intensity);
  
  //apply a threshold to make the simulation a bit more realistic
  for(int i=0; i<n_x; i++){
    for(int j=0; j<n_y; j++){
      intensity.set(i,j,intensity.get(i,j)-noise_level);
      if(intensity.get(i,j)<0)
        intensity.set(i,j,0);
    }
  }

  //write the output to file (use log scale)
  write_tiff("part_char_real_sim_intensity.tiff",intensity,false);
  write_tiff("part_char_log_sim_intensity.tiff",intensity,true);
  
  write_dbin("image_files/part_char_sim.dbin",intensity);

  /******** get the support from file ****************************/

  Double_2D support(n_x,n_y);
  status = read_tiff(support_file_name, support);

  /*************** do the reconstruction *******************/

  //create a project object and set the options.
  my_partial.set_support(support);
  my_partial.set_intensity(intensity);
  my_partial.set_algorithm(HIO);

  //set the inital guess to be random inside the support
  //and zero outside. Note that this must be called
  //after "my_partial.set_support()"
  my_partial.initialise_estimate(0);
  
  //make a temporary arrary
  Double_2D result(n_x,n_y);

  //apply the projection operators
  for(int i=0; i<hio_iterations; i++){

    cout << "iteration " << i << endl;

    my_partial.iterate();

    if(i%output_iterations==0){

      ostringstream temp_str ( ostringstream::out ) ;
      first_guess.get_2d(MAG,result);
      temp_str << "sim_result_" << i << ".ppm";
      write_ppm(temp_str.str(), result);

      my_partial.apply_shrinkwrap();

    }
  }

  my_partial.set_algorithm(ER);

  for(int i=hio_iterations; i<(er_iterations+hio_iterations+1); i++){
    
    cout << "iteration " << i << endl;
    
    my_partial.iterate();

    if(i%output_iterations==0){

      ostringstream temp_str ( ostringstream::out ) ;
      first_guess.get_2d(MAG,result);
      temp_str << "sim_result_" << i << ".ppm";
      write_ppm(temp_str.str(), result);

      my_partial.apply_shrinkwrap();

    }
  }



  
  return 0;
}

