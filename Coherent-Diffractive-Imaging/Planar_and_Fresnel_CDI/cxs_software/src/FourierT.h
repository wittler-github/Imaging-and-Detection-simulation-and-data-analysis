/**
 * @file FourierT.h
 * @class FourierT
 * @author Nadia Davidson
 * @date Last modified on 6/1/2011
 *
 * @brief An interface to the fftw3 algorithms
 * 
 * This class is now only provided for the user's own need. All
 * ffts in the other classes of this library are preformed
 * within the Complex_2D class, removing the need to copy 
 * large arrays around.
 *
 * This class provides an interface to the fftw3 algorithms.  This
 * makes the reconstruction code a nit neater as the creation and
 * destruction of fftw plans is kept track of here. It also provides
 * methods for copying between Complex_2D objects and the fftw arrays
 * used by the fourier transforming algorithm.
 */

#ifndef FOURIERT_H
#define FOURIERT_H

#include <fftw3.h>

#define FAILURE 0
#define SUCCESS 1

#define REAL 0
#define IMAG 1
#define MAG 2
#define PHASE 3

class Complex_2D;

class FourierT{

 private:
  
  /* The size of the array in x */
  int nx;

  /* The size of the array in y */
  int ny;

  /* A complex array which holds the input to the transform */
  fftw_complex * original;

  /* A complex array which holds the output from the transform */
  fftw_complex * transformed;

  /* A fftw plan for forward fourier transforms */
  fftw_plan f_forward;

  /* A fftw plan for backward fourier transforms */
  fftw_plan f_backward;

 public:

  /**
   * Create a FourierT for Complex_2Ds with dimensions of x_size and
   * y_size. An instance of this class maybe reused for any Complex_2D
   * with these dimensions. The fftw plans are set up with the
   * "MEASURED" option.
   *
   * @param x_size The size of the array in x.
   * @param y_size The size of the array in y.
   */
  FourierT(int x_size, int y_size);

  /**
   * Destructor
   */
  ~FourierT();

  /**
   * Forward fourier transform the given Complex_2D object. The
   * Complex_2D is scaled to give the same normalisation as before the
   * transform.
   *
   * @param c_in The Complex_2D to transform. The object will be
   * modified with the result.
   */
  void perform_forward_fft(Complex_2D & c_in);

  /**
   * Backward fourier transform the given Complex_2D object. The
   * Complex_2D is scaled to give the same normalisation as before the
   * transform.
   *
   * @param c_in The Complex_2D to transform. The object will be
   * modified with the result.
   */
  void perform_backward_fft(Complex_2D & c_in);
  
 private:

  /**
   * Used by the perform_forward_fft and perform_backward_fft
   * functions. The fast memcpy method is used.
   */
  void copy_to_fftw_array(fftw_complex * array , Complex_2D & c);

  /**
   * Used by the perform_forward_fft and perform_backward_fft
   * functions. The array is scaled while copying to ensure that the
   * result has the same normalisation as it did prior to the
   * transform.
   */
  void copy_from_fftw_array(fftw_complex * array , Complex_2D & c);
  
};

#endif
