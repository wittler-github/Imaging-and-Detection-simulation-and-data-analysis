// Copyright 2011 Nadia Davidson 
// for The ARC Centre of Excellence in Coherent X-ray Science. 
//
// This program is distributed under the GNU General Public License. 
// We also ask that you cite this software in publications where you made 
// use of it for any part of the data analysis.
//
// date last modified: 21/06/2013

/**
 * @file PlanarCDI_simulation_example.c
 *
 * \a workbook template for PlanarCDI_simulation_example.c.
 *
 */

#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib> 
#include <io.h>
#include <Complex_2D.h>
#include <Double_2D.h>
#include <PlanarCDI.h>
#include <sstream>

using namespace std;


/**************************************/
int main(void){

/****** get the object from an image file ****************/

  int n_x=1024;
  int n_y=1024;

//get the data from file
  Double_2D object;

  read_tiff("object_file_name", object);

  Complex_2D object_estimate(n_x, n_y);
  Complex_2D input(n_x,n_y);

  for(int i=0; i < n_x; i++){
    for(int j=0; j < n_y; j++){
      input.set_value(i, j, REAL, 1.0/sqrt(2.0)*object.get(i,j));
      input.set_value(i, j, IMAG, 1.0/sqrt(2.0)*object.get(i,j));
    }
  }

  PlanarCDI planar(object_estimate);

  planar.propagate_to_detector(input);

  Double_2D object_mag(n_x,n_y);
  
  input.get_2d(MAG_SQ,object_mag);
  write_tiff("magnitude_file_name.tif", object_mag );
 
  return 0;
}

