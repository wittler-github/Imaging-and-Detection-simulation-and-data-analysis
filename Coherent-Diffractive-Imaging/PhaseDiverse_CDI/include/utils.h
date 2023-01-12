// Copyright 2011-2013 Nadia Davidson & Daniel Rodgers-Pryor for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <Double_2D.h>

#define DEFAULT_GAUSSIAN_KERNEL_SIZE_IN_STD_DEVIATIONS 3 // Use a matrix that is (at least) this many std deviations wide/high to discretly represent gaussians


//class Double_2D;
class Complex_2D;

void crop(Double_2D & image, Double_2D & new_image, int x_start, int y_start);
void rescale(Double_2D & image, double scale);
void slow_align(Double_2D & image1, Double_2D & image2, int & offset_x, int & offset_y, int step_size=8.0, int min_x=0, int max_x=0, int min_y=0, int max_y=0);

void align(Double_2D & first_image, Double_2D & second_image,
	int & offset_x, int & offset_y,
  int min_x=0, int max_x=0,
	int min_y=0, int max_y=0,
	Double_2D * first_image_weights = 0,
	Double_2D * second_image_weights = 0,
	double overlap_fraction = 0.2);

double edges(Double_2D & image);
double line_out(Double_2D & image);
double calculate_high_frequency_ratio(Double_2D & image);
double calculate_image_entropy(Double_2D & image);
double calculate_gradients(Double_2D & image, double threshold=0.03);
double calculate_mean_difference(Double_2D & image);
double calculate_average_energy_density(Double_2D & image);
double vollaths_4(Double_2D & image);
double vollaths_5(Double_2D & image);
double deviation_from_zero(Double_2D & image);
double count_pixels(Double_2D & image, double threshold);
double diff_of_squares(Double_2D & image1, Double_2D & image2);
double simple(Double_2D & image, double scale);
double sobel_gradient(Double_2D & image);
double laplace_gradient(Double_2D & image);
double edge_grad(Double_2D & image, Double_2D & mask);
double calculate_image_entropy_2(Double_2D & image);

void interpolate( const Complex_2D & original, Complex_2D & big);
void interpolate( const Double_2D & original, Double_2D & big);
void shrink( const Complex_2D & original, Complex_2D & small);
void shrink( const Double_2D & original, Double_2D & small);

//Legendre related functions
Double_2D  legroots(double n);
Double_2D fill_legmatrix(std::vector<double>, int norder);

void solve_gep(Complex_2D & A, Complex_2D & B, std::vector<double> & eigen);

void* smalloc(size_t size);
double sq(double x);
double avg(double a, double b);
int fuzzy_eq(double a, double b, double tolerance);
double* get_gaussian_vector(double std_dev, unsigned int n);
Double_2D vector_convolution(Double_2D const & m, double* v, unsigned int n, bool rowwise);
Double_2D gaussian_convolution(Double_2D const & m, double lx, double ly,
  double kernel_x_size_in_std_dev=DEFAULT_GAUSSIAN_KERNEL_SIZE_IN_STD_DEVIATIONS, double kernel_y_size_in_std_dev=DEFAULT_GAUSSIAN_KERNEL_SIZE_IN_STD_DEVIATIONS);
Double_2D radial_gaussian_convolution(Double_2D const & m, double lr, double kernel_size_in_std_dev=DEFAULT_GAUSSIAN_KERNEL_SIZE_IN_STD_DEVIATIONS);
void convolve(Double_2D & array, double gauss_width, int pixel_cut_off);

/**
 * Template class for a generic function of the form:
 *   f:R->R
 *
 * An instance of a subclass of this must be used as the argument to minimise_function
 */
class MathFunction {
  public:
      virtual double call(double x) = 0;
};

// Find minima of f, within bounds, to within tolerance
double minimise_function(MathFunction & f, double left, double guess, double right, double tolerance);

#endif

