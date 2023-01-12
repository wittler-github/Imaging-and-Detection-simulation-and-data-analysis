// Copyright 2011 Nadia Davidson 
// for The ARC Centre of Excellence in Coherent X-ray Science. 
//
// This program is distributed under the GNU General Public License. 
// We also ask that you cite this software in publications where you made 
// use of it for any part of the data analysis.

/**
 * @file PlanarCDI_example.c
 *
 * \a PlanarCDI_example.c This example reconstructs some planar
 * diffraction data (Lachie's data). The shrinkwrap algorithm is used
 * to improve the reconstruction. A combination of HIO and the
 * error-reduction algorithm are used.
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include <io.h>
#include <Complex_2D.h>
#include <PlanarCDI.h>
#include <Double_2D.h>
#include <iomanip>

using namespace std;

int main(void){

  int nx = 1024;
  int ny = 1024;

  int n_iterations = 100;

  Double_2D data, support;
  read_image("magnitude_file_name.tif", data);
  read_image("support_file_name.tif", support);

  Complex_2D object_estimate(nx, ny);
  PlanarCDI planar(object_estimate);

  planar.set_support(support);
  planar.set_intensity(data);
  planar.initialise_estimate(0);

  for(int i=0; i <n_iterations; i++){
    planar.iterate();

  }

  write_cplx("result.cplx", object_estimate);
  Double_2D object_mag(nx, ny);
  object_estimate.get_2d(MAG, object_mag);
  write_image("magnitude_data_file.tiff", object_mag);

  return 0;
}

