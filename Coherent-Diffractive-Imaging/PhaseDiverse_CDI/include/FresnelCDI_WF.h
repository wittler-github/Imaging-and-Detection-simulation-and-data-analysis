/**
 * @file FresnelCDI_WF.h
 * @class FresnelCDI_WF
 * @author  Nadia Davidson <nadiamd@unimelb.edu.au> 
 * @data Last modified on 7/1/2011
 *
 * @brief The class which performs Fresnel CDI white field
 * reconstruction.
 *
 * This is the class used for reconstructing white-field data prior to
 * running Fresnel reconstruction with FresnelCDI. It inherits most
 * methods from BaseCDI, so please look at the documentation of this
 * class also. Although there are some differences in the underlying
 * code between this class and the planar case, the interface is
 * generally unchanged. Therefore users should refer to the
 * instructions for BaseCDI to understand how to use a FresnelCDI_WF
 * object in their own code. Only the differences relevant to users
 * will be documented here.
 *
 */

#ifndef FCDI_WF_H
#define FCDI_WF_H

#include "BaseCDI.h"


//forward declarations
class Complex_2D;

class FresnelCDI_WF: public BaseCDI{

 protected:

  /** the beam wavelength */
  double wavelength;

  /** the distance from the zone plate to the focal point */
  double zone_to_focal_length; 

  /** the distance from the focal point to the detector */
  double focal_to_detector_length;

  /** the pixel size */
  double pixel_length;

   /** an array which holds a constants we use when propagating between
       difference planes */
  Complex_2D coefficient;
  
  /** an array which holds a constants we use when propagating between
      difference planes */
  //  Complex_2D backward_coefficient;
  
 public:
  
  /** 
   * Create a FresnelCDI_WF object. All of the dimensions which are given as
   * input should be in the same units. e.g metres.
   *
   * @param initial_guess The complex 2-D field which is modified by
   * this object. This represents the illuminating white-field in the
   * detector plane. It may be pre-initialised to a best first guess
   * (e.g manually); it may be loaded from a file (e.g. the output
   * from a previous reconstruction job); or it may be unitialised,
   * and the default initialisation provided by FresnelCDI_WF may be
   * used.
   * @param beam_wavelength Wavelength of the beam 
   * @param zone_focal_length The distance from the zone plate to the
   * focal point
   * @param focal_detector_length The distance from the focal point to the
   * detector
   * @param pixel_size The diameter of one of the detector pixels.
   */
  FresnelCDI_WF(Complex_2D & initial_guess,
		double beam_wavelength,
		double zone_focal_length,
		double focal_detector_length,
		double pixel_size,
		int n_best=1);
  
  
  /**
   * The destructor for this class
   */
  virtual ~FresnelCDI_WF(){};


  /**
   * The iterate() method of BaseCDI is overridden here with 3-plane
   * propagation.
   */
  int iterate();  

  /**
   * Propagate to the zone plate plane using the paraxial 
   * free-space equation.
   * @param c The Complex_2D field to propagate
   */
  virtual void propagate_from_detector(Complex_2D & c);

  /**
   * Propagate to the detector plane using the paraxial free-space
   * equation.
   * @param c The Complex_2D field to propagate
   */
  virtual void propagate_to_detector(Complex_2D & c);


  /**
   * This method overrides the one in BaseCDI. The white-field is
   * initialised to be zero outside the support (note that we assume
   * here that the support for the zone plate is approx. the same size
   * and shape as the support region of the white-field in the detector
   * plane. Inside the support region, the real component is set to
   * the random number and the imaginary component is set to zero
   * (this makes the reconstruction much fast than it would be if both
   * components were random).
   *
   * @param seed The seed for the random number generator
   * 
   */
  virtual void initialise_estimate(int seed=0);

  /**
   * This method is an alternative to setting the support using an
   * image. You should pass the zone-plate diameter and (optionally)
   * how much bigger the support should be than the zone-plate. The
   * support will be set to a circle, centered at the center of the
   * zone-plate plane, with a diameter equal to the zone plate
   * diameter time the size parameter. The support edge is softened by
   * convoluting it with a gaussian function (width=3 pixels).
   *
   * @param z_diameter Diameter of the zone-plate (given in the same
   * units as the dimensions given to the constructed).  

   * @param size What proportion of the zone-plate plane will be used
   * for the support, by default this is 1.01 (the support is 1%
   * larger than the support in diameter).
   */
  virtual void set_support(double z_diameter, double size=1.01);


  /**
   * This method is not available for FresnelCDI_WF.
   */
  void set_algorithm(int alg){
    std::cout << "WARNING: Algorithm can not be set when performing Fresnel "
	      << "CDI white field recovery" << std::endl;
  }
  
  /**
   * This method is not available for FresnelCDI_WF.
   */
  void print_algorithm(){
    std::cout << "Using the default Fresnel CDI "
	      << "white field recovery algorithm" << std::endl;
  }

  /**
   * This method is not available for FresnelCDI_WF.
   *
   */
  void set_relaxation_parameter(double relaxation_parameter){
    std::cout << "WARNING: The relaxation parameter can not "
	      << "set when performing Fresnel white-field "
	      << "reconstruction" << std::endl;

  };

  /**
   * This method is not available for FresnelCDI_WF.
   *
   */
  void set_custom_algorithm(double m1, double m2, double m3, double m4, 
			    double m5, double m6, double m7, double m8, 
			    double m9, double m10){
  };

  
  void multiply_factors(Complex_2D & c, int direction);

};

#endif
