// Copyright 2013 Daniel Rodgers-Pryor for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

#include <cmath>
#include <stdlib.h>
#include <Double_2D.h>
#include <Complex_2D.h>
#include <sstream>
#include <utils.h>
#include <PartialCharCDI.h>
#include <cstdlib> 

using namespace std;


PartialCharCDI::PartialCharCDI(Complex_2D & initial_guess, 
			      double z_dist, 
			      double beam_energy, 
			      double pixel_x_size, 
			      double pixel_y_size, 
			      unsigned int n_best)
  :BaseCDI(initial_guess,n_best),
    iteration(0),
    lx(DEFAULT_INITIAL_LY),
    ly(DEFAULT_INITIAL_LX),
    wavelength(plank_h*speed_c/beam_energy),
    z(z_dist),
    px_x_size(pixel_x_size),
    px_y_size(pixel_y_size),
    measured_intensity(),
    minima_search_bounds_coefficient(DEFAULT_MINIMA_SEARCH_BOUNDS_COEFFICIENT),
    minima_search_tolerance(DEFAULT_MINIMA_SEARCH_TOLERANCE),
    minima_recalculation_interval(DEFAULT_MINIMA_RECALCULATION_INTERVAL),
    minima_moving_average_weight(DEFAULT_MINIMA_MOVING_AVERAGE_WEIGHT){}


/**
 * Default: 0.7, 0.7
 *
 * Set an initial guess for the beam coherence lengths (in pixels).
 *
 * Note that while bad values for this guess won't stop convergence, they can slow it down significantly.
 */
void PartialCharCDI::set_initial_coherence_guess(double lx_, double ly_){
  // Invert to FT into fourier space
  lx = 1.0/lx_;
  ly = 1.0/ly_;
}

/**
 * Set an initial guess for the beam coherence lengths in m, measured in the object plane.
 */
void PartialCharCDI::set_initial_coherence_guess_in_m(double lx_, double ly_){
  lx = 2*px_x_size/(lx_*z*wavelength);
  ly = 2*px_y_size/(ly_*z*wavelength);
}

/**
 * Default = 2.5
 *
 * Set the factor that will be used to determine the search region.
 * Search will go from 1/coef to coef times the current estimate.
 *
 * If you are very unsure of your initial guess of lx and ly (beam coherence lengths) 
 * or if you are just using the default guess of 0.7 pixels in each direction, then you might want to
 * increase this to 3.
 */
void PartialCharCDI::set_minima_search_bounds_coefficient(double coef){
  minima_search_bounds_coefficient = coef;
}

/**
 * Default = 0.1
 *
 * The the tolerance of the lx/ly minima search in pixels (search will terminate when uncertainty is less than this value).
 */
void PartialCharCDI::set_minima_search_tolerance(double tol){
  // Invert to FT into fourier space
  minima_search_tolerance = 1.0/tol;
}

/**
 * The the tolerance of the lx/ly minima search in object-plane-meters (search will terminate when uncertainty is less than this value).
 */
void PartialCharCDI::set_minima_search_tolerance_in_m(double tol){
  double px_radial_size = avg(px_x_size, px_y_size); // Use an average pixel size since this tolerance will apply to both dimensions
  minima_search_tolerance = 2*px_radial_size/(tol*z*wavelength);
}

/**
 * Default = 0.01
 *
 * lx and ly values will be time averaged with the previous estimate (to smooth
 * the slightly-erratic changes from iteration to iteration).
 * w is the weight of the newly calculated value. The current estimate will be weighted at (1-w) to
 * form a new current estimate.
 * This is an exponential moving average.
 *
 * Setting w to 1.0 will be the same as having no moving average.
 */
void PartialCharCDI::set_minima_moving_average_weight(double w){
  minima_moving_average_weight = min(max(w, 0.0), 1.0); // Ensure that the weight is in the range [0,1]
}

/** 
 * Default = 5
 *
 * The values of lx and ly will only be updated every ival iterations (since the process is rather slow).
 * There's not much advantage to doing this more frquently than the default of 5 unless there are very few
 * total iterations being calculated.
 */
void PartialCharCDI::set_minima_recalculation_interval(unsigned int ival){
  minima_recalculation_interval = ival;
}

/**
 * Retrieve the calculated coherence length of the beam (the std. deviation of the fitted gaussian) in 
 * the x direction, in meters, in the object plane.
 *
 * Call this after iterating enough to produce a satisfactory image.
 * This value will only be accurate to within the tolerances set by set_minima_search_tolerance_in_m.
 */
double PartialCharCDI::get_x_coherence_length(){
  return 0.5*z*wavelength/(lx*px_x_size);
}

/**
 * Retrieve the calculated coherence length of the beam (the std. deviation of the fitted gaussian) in 
 * the y direction, in meters, in the object plane.
 *
 * Call this after iterating enough to produce a satisfactory image.
 * This value will only be accurate to within the tolerances set by set_minima_search_tolerance_in_m.
 */
double PartialCharCDI::get_y_coherence_length(){
  return 0.5*z*wavelength/(ly*px_y_size);
}

/**
 * Retrieve the calculated x coherence length of the beam in pixels. This can be useful for interacting
 * simulated data on a pixel-scale.
 */
double PartialCharCDI::get_x_coherence_length_in_pixels(){
  // Invert to FFT back into real space
  return 1.0/lx;
}

/**
 * Retrieve the calculated y coherence length of the beam in pixels. This can be useful for interacting
 * simulated data on a pixel-scale.
 */
double PartialCharCDI::get_y_coherence_length_in_pixels(){
  // Invert to FFT back into real space
  return 1.0/ly;
}


void PartialCharCDI::initialise_estimate(int seed){
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

void PartialCharCDI::set_intensity(const Double_2D &detector_intensity){
  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){
      intensity_sqrt.set(i,j,sqrt(detector_intensity.get(i,j)));
    }
  }

  // Store a copy of the un-square-rooted intensity values:
  measured_intensity = detector_intensity;
}

void PartialCharCDI::propagate_to_detector(Complex_2D & c){
    c.perform_forward_fft();
    c.invert(true);
}

void PartialCharCDI::propagate_from_detector(Complex_2D & c){
  c.invert(true);
  c.perform_backward_fft(); 
}


//this overwrites the function of the same
//name in BaseCDI
void PartialCharCDI::scale_intensity(Complex_2D & c){
  double norm2_mag=0;
  double norm2_diff=0;
  double current_int_sqrt=0;
  double current_mag=0;

  Double_2D estimate_intensity(nx, ny);
  c.get_2d(MAG_SQ, estimate_intensity); // Get the intensity of the current estimate
  Double_2D optimally_convoluted_estimate;
  
  if(iteration % minima_recalculation_interval == 0){ // Every few iterations:
    // Recalculate the optimal gaussian convolution of the estimated intensity to minimise the difference with the measured intensity:
    optimally_convoluted_estimate = get_convoluted_intensity_estimate(estimate_intensity, measured_intensity); // lx and ly attributes are updated in here
  } else{ // Use last lx and ly values:
    optimally_convoluted_estimate = gaussian_convolution(estimate_intensity, lx, ly);
  }
  optimally_convoluted_estimate.sq_root(); // get sqrt

  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){
      // Reset the magnitude to match de-convoluted measurement intensities:
      if(beam_stop==0 || beam_stop->get(i,j)>0){
        current_mag = c.get_mag(i, j);

        current_int_sqrt = current_mag * intensity_sqrt.get(i,j) / optimally_convoluted_estimate.get(i,j); // Softened intensity-constraint
        c.set_mag(i,j,current_int_sqrt);
        
        //calculate the error
        norm2_mag += current_int_sqrt*current_int_sqrt;
        norm2_diff += (current_mag-current_int_sqrt)
          *(current_mag-current_int_sqrt);
      }
    }
  }
  current_error = (norm2_diff/norm2_mag);

  iteration++; // Increment iteration count
}


/**
 * Find a gaussian matrix g with standard deviations in x and y that minimise
 * the energy difference between the gaussian-convoluted current estimated 
 * intensity and the measured intensity.
 * Returns the convolution of c and g (ie. the gaussian blurred estimate that best
 * matches the measured intensity). This can be thought of as a coherence-corrected
 * estimate.
 *
 * See: Applied Physics Letters 99 (2011) - 'Simultaneous sample and spatial
 * coherence characterisation using diffractive imaging' [doi: 10.1063/1.3650265]
 * for further details.
 */
Double_2D PartialCharCDI::get_convoluted_intensity_estimate(Double_2D const & estimated_intensity, Double_2D const & measured_intensity){
  double lbound, rbound;

  // Find optimal lx estimate with exponential moving average:
  // Use last values as initial guess and a fixed multiple/division of the last value as the search bounds
  lbound = max(MIN_SEARCH_LBOUND, lx / minima_search_bounds_coefficient); // Don't search the region where the gaussian is delta-function like
  rbound = max(MIN_SEARCH_RBOUND, lx * minima_search_bounds_coefficient); // Dont let the bound size get too small
  error_in_lx f_lx(measured_intensity, estimated_intensity, ly); // Error as a function of lx (using the previous value of ly w.l.o.g.)
  
  // Find new error-minimising value and time-average with old value:
  lx = ((1-minima_moving_average_weight) * lx) + (minima_moving_average_weight * minimise_function(f_lx, lbound, lx, rbound, minima_search_tolerance));

  // Find optimal ly estimate with exponential moving average:
  lbound = max(MIN_SEARCH_LBOUND, ly / minima_search_bounds_coefficient);
  rbound = max(MIN_SEARCH_RBOUND, ly * minima_search_bounds_coefficient);
  error_in_ly f_ly(measured_intensity, estimated_intensity, lx); // Error as a function of ly

  // Find new error-minimising value and time-average with old value:
  ly = ((1-minima_moving_average_weight) * ly) + (minima_moving_average_weight * minimise_function(f_ly, lbound, ly, rbound, minima_search_tolerance));
  
  return gaussian_convolution(estimated_intensity, lx, ly); // Return the convolution of estimated_intensity with a gaussian of optimal parameters
}

/**
 * Returns the energy difference between the measureed intensity, and the intensity
 * estimate assuming partial coherence described by gaussian std deviation parameters lx & ly
 */
double PartialCharCDI::convoluted_estimate_error(Double_2D const & measured_intensity, Double_2D const & estimated_intensity, double lx, double ly){
  Double_2D convoluted_estimate = gaussian_convolution(estimated_intensity, lx, ly);

  // Subtract the convoluted_estimate from measured_intensity in preparation for summing differences:
  convoluted_estimate.scale(-1);
  convoluted_estimate.add(measured_intensity);
  
  // Return the sum of absolute differences between convoluted_estimate and measured_intensity:
  return convoluted_estimate.get_abs_sum();
}
