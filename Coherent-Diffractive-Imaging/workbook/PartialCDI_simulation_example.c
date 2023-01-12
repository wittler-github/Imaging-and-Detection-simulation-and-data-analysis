// Copyright 2012 T'Mir Juliusn 
// for The ARC Centre of Excellence in Coherent X-ray Science. 
//
// This program is distributed under the GNU General Public License. 
// We also ask that you cite this software in publications where you made 
// use of it for any part of the data analysis.
//
// date last modified: 05/07/2013

/**
 * @file PartialCDI_simulation_example.c
 *
 * \a PartialCDI_simulation_example.c template. The object from the
 * PartialCDI_example (in fact I used the reconstructed image as the
 * object) is used to simulate a diffraction pattern for partial
 * CDI. The diffraction pattern is thresholded to make it more
 * realistic and then CDI reconstruction is performed.
 *
 */

#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib> 
#include <io.h>
#include <Complex_2D.h>
#include <Double_2D.h>
#include <PartialCDI.h>
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

  //The number of legendre polynomials must be greater than or equal to the 
  //number of modes
  int nleg=8;
  int nmodes=7;


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

  //fill the complex no. with image data
  Complex_2D input(n_x,n_y);
  for(int i=0; i<n_x; i++){
    for(int j=0; j<n_y; j++){
      input.set_value(i,j,REAL, 1.0/sqrt(2.0)*(data.get(i,j)));
      input.set_value(i,j,IMAG, 1.0/sqrt(2.0)*(data.get(i,j)));
    }
  }

  /**** create the projection/reconstruction object *****/

  //Coherence lengths x and y in m
  double lcx = 13.3e-6;
  double lcy = 40.0e-3;

  //Pixel size detector in m
  double psize_x=13.5e-6;
  double psize_y=13.5e-6;

  //Energy of the beam in eV
  double e_beam=1400.0;

  //Distance between detector and sampl in metres
  double z_sd=1.4;

  PartialCDI my_partial(input, lcx, lcy, psize_x, psize_y, e_beam, z_sd, nleg, nmodes);

  Complex_2D pattern(n_x,n_y);

  my_partial.set_threshold(+1.0e-6); //Set the threshold as a fraction of the dominant mode.

  Double_2D result(n_x, n_y);

  //Uncomment the section below to output images of each of the included modes
  /*for(int i=0; i< nmodes*nmodes; i++){
    ostringstream temp_str0 ( ostringstream::out ) ;
    my_partial.get_mode(i).get_2d(REAL,result);
    temp_str0 << "mode"<<i<<".ppm";
    write_image(temp_str0.str(), result);
  }
  */

  //propagate to the detector plane
  my_partial.set_transmission(input);

  Double_2D intensity=my_partial.propagate_modes_to_detector();

  //apply a threshold to make the simulation a bit more realistic
  for(int i=0; i<n_x; i++){
    for(int j=0; j<n_y; j++){
      intensity.set(i,j,intensity.get(i,j)-noise_level);
      if(intensity.get(i,j)<0)
	intensity.set(i,j,0);
    }
  }

  return 0;
}

