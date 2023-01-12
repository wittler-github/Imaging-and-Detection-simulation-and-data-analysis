// Copyright 2013 Daniel Rodgers-Pryor for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

#ifndef PCCDI_H
#define PCCDI_H

#define _USE_MATH_DEFINES

#include "BaseCDI.h"
#include "utils.h"
#include <map>

// Initial guess for coherence-spread (measured in pixels)
#define DEFAULT_INITIAL_LX 0.7
#define DEFAULT_INITIAL_LY 0.7

#define DEFAULT_MINIMA_SEARCH_BOUNDS_COEFFICIENT 2.5 // Use ([current_estimate]/[this], [current_estimate] * [this]) as the search interval for the energy-minimising beam coherence value
#define DEFAULT_MINIMA_SEARCH_TOLERANCE 0.1 // Tolerated error in finding the energy-minimising coherence values
#define DEFAULT_MINIMA_MOVING_AVERAGE_WEIGHT 0.01 // The minima will be updated to be: this*new_measurement + (1-this)*last_value
#define DEFAULT_MINIMA_RECALCULATION_INTERVAL 5 // The minima will be updated once every [this many] generations

#define MIN_SEARCH_LBOUND 0.2 // Don't bother searching near zero where the gaussian is delta-fn-like
#define MIN_SEARCH_RBOUND 1.0 // Don't let the search region get too small
    
const double plank_h = 4.13566733e-15;
const double speed_c = 2.99792458e8;


class PartialCharCDI:public BaseCDI {
  public:
    PartialCharCDI(Complex_2D & initial_guess, double z_dist=2, double beam_energy=1.239841e-6, double pixel_x_size=1, double pixel_y_size=1, unsigned int n_best=0);
    
    // Accessor functions to alter initial parameters:
    void set_initial_coherence_guess(double lx_, double ly_);
    void set_initial_coherence_guess_in_m(double lx_, double ly_);
    virtual void initialise_estimate(int seed=0);
    void set_intensity(const Double_2D &detector_intensity);
    
    // Accessor functions to alter search parameters:
    void set_minima_search_bounds_coefficient(double coef);
    void set_minima_search_tolerance(double tol);
    void set_minima_search_tolerance_in_m(double tol);
    void set_minima_moving_average_weight(double w);
    void set_minima_recalculation_interval(unsigned int ival);

    // Result Accessor functions:
    double get_x_coherence_length();
    double get_y_coherence_length();
    double get_x_coherence_length_in_pixels();
    double get_y_coherence_length_in_pixels();

    virtual void propagate_to_detector(Complex_2D & c);
    virtual void propagate_from_detector(Complex_2D & c);
    
    void scale_intensity(Complex_2D & c); // Scale the magnitude of the complex estimate using a convolution of estimate intensity

    // Allow these MathFunctions to access convoluted_estimate_error externally:
    friend class error_in_lx;
    friend class error_in_ly;

  protected:
    /**
     * Last calculated coherence lengths.
     *
     * Note: 'coherence length' is defined as the standard deviation of the coherence function in the DETECTOR-PLANE. This means
     * that these are the fourier transform (ie. inverse since it's gaussian) of the coherence length in the object-plane.
     * */
    double lx;
    double ly;

    double z; // Distance from sample to detector
    double wavelength; // Wavelength of the beam

    // Size of detector pixels in meters
    double px_x_size;
    double px_y_size;

    Double_2D measured_intensity; // Un-square-rooted detector intensity - needed for measuring gaussian convolution error
    
  private:
    unsigned int iteration; // Count of the number of iterations in the scale_intensity function (to know when to reevaluate lx and ly)
    unsigned int minima_recalculation_interval; // Recalculate lx and ly once every minima_recalculation_interval iterations
    
    double minima_search_bounds_coefficient; // Search bounds for lx, ly minima will be (current_value/coefficient, current_value*coefficient)
    double minima_search_tolerance; // Precision to which lx and ly are calculated
    double minima_moving_average_weight; // A running average will be used to estimate lx and ly. Each new update will be weighted at this fraction.
    
    Double_2D get_convoluted_intensity_estimate(Double_2D const & estimated_intensity, Double_2D const & measured_intensity);
    
    /**
     * Take the integral over the absolute difference between the measured_intensity and the convolution of estimated_intesity
     * and a gaussian of standard deviations lx and ly in x and y directions respectivley.
     *
     * This is a measure for how accurate the estimate is given known coherence lengths lx and ly. It is used internally to
     * automatically fit for lx and ly (by minimising its value in lx and ly) every few iterations.
     *
     * This call is quite expensive.
     */
    static double convoluted_estimate_error(Double_2D const & measured_intensity, Double_2D const & estimated_intensity, double lx, double ly);
};


/*
 * Tools for finding minima of lx and ly:
 */

class error_in_lx : public MathFunction {
  // convoluted_estimate_error as a function only of lx (with caching)
  public:
      error_in_lx(Double_2D const & measured_intensity, Double_2D const & estimated_intensity, double ly) : estimated_intensity_(estimated_intensity), measured_intensity_(measured_intensity), ly_(ly) {
        cache.clear();
      }
      virtual double call(double lx){
        double result;
        if(cache.find(lx) == cache.end()){ // Cache the expensive call
          result = PartialCharCDI::convoluted_estimate_error(measured_intensity_, estimated_intensity_, lx, ly_);
          cache[lx] = result;
        } else{
          result = cache[lx];
        }

        return result;
      }
  private:
          double ly_;
          Double_2D const & estimated_intensity_;
          Double_2D const & measured_intensity_;
          std::map<double, double> cache;
};

class error_in_ly : public MathFunction {
  // convoluted_estimate_error as a function only of ly (with caching)
  public:
      error_in_ly(Double_2D const & measured_intensity, Double_2D const & estimated_intensity, double lx) : estimated_intensity_(estimated_intensity), measured_intensity_(measured_intensity), lx_(lx) {
        cache.clear();
      }
      virtual double call(double ly){
        double result;
        if(cache.find(ly) == cache.end()){ // Cache the expensive call
          result = PartialCharCDI::convoluted_estimate_error(measured_intensity_, estimated_intensity_, lx_, ly);
          cache[ly] = result;
        } else{
          result = cache[ly];
        }

        return result;
      }
  private:
          double lx_;
          Double_2D const & estimated_intensity_;
          Double_2D const & measured_intensity_;
          std::map<double, double> cache;
};

#endif
