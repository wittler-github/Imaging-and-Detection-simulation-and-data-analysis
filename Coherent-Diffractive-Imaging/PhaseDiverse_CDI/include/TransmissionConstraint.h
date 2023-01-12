// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

/**
 * @file TransmissionConstraint.h
 * @author  Nadia Davidson <nadiamd@unimelb.edu.au> 
 */


#ifndef TRANS_H
#define TRANS_H
#include <cstdlib>
#include <vector>
#include <iostream>
#include "Double_2D.h"

//class Double_2D;
class Complex_2D;


/**
 * @class ComplexConstraint
 * @author  Nadia Davidson <nadiamd@unimelb.edu.au> 
 *
 * @brief Create a constraint on the refractive index in
 *        one region of the reconstruction.
 *
 * This class is based on the paper "Use of a complex constraint in
 * coherent diffractive imaging", J. N. Clark et. al., 2010. The
 * notation used there is also used here. Users should have knowledge
 * of this or similar papers before using this class.
 * 
 * A section of the reconstructed transmission function will be
 * updated according to a restriction on the value of c (c =
 * beta/delta).  If the material is of a known element, then c can be
 * fixed to the know value, and the phase and magnitude of the
 * transmission function updated accordingly. Alternatively, c can be
 * left to float, and calculated for each iteration from a mean over
 * the defined region. The parameters alpha1 and alpha2 are used to
 * control the strength of the constraint.
 * 
 * Instances of this class should be pass to a TransmissionConstraint
 * object. If a sample consists of different regions which are locally
 * homogeneous in element, multiple instances of this class should be
 * created, each one defining a different region of the sample, and
 * then each should be added to the TransmissionConstraint object.
 */
class ComplexConstraint{

 protected:
  
  /** constraint strength parameter for the amplitude */
  double alpha1;

  /** constraint strength parameter for the phase */
  double alpha2;

  /** region over which the constraint should be applied */
  Double_2D * region;

  /** the value of c = beta / delta */
  double c_mean;

  /** flag as to whether c is fixed or allowed to be updated */
  bool fixed_c;

 public:

  /**
   * Construct a complex constraint.
   *
   * @param region A Double_2D array with positive element values
   * defining the region for which this constraint should be applied.
   * Once set, it can not be reset.
   * @param alpha1 The constraint strength parameter for the magnitude.
   * @param alpha2 The constraint strength parameter for the phase.
   */
  ComplexConstraint(Double_2D & region,
		   double alpha1,
		   double alpha2){
    this->alpha1 = alpha1;
    this->alpha2 = alpha2;
    this->region = &region;
    fixed_c = false;
  };

  /**
   * Set a fixed value of c = beta/delta. If set, you will no longer
   * be able to use the set_c_mean() method. This value will be used
   * for all iterations.
   *
   * @param c_mean The value of c = beta/delta.
   */
  void set_fixed_c(double c_mean){
    fixed_c = true;
    this->c_mean = c_mean;
  };

  /**
   * Set the value C = beta/delta. This function is used by
   * TransmissionConstraint to set the value of C, and does probably
   * not need to be called by a user. The function
   * TransmissionConstraint::apply_constraint() calculated C as:
   *
   * \f[ 
   * \bar{C} = \frac{\sum{lnA}}{\sum{\phi}}
   * \f]
   *
   * Where A is the transmission function amplitude, \f$\phi\f$ is
   * the phase and the summation runs over all values in the
   * transmission function array for the defined region.
   *
   * @param c_mean The value of C = beta/delta.
   */
  void set_c_mean(double c_mean){
    if(!fixed_c)
      this->c_mean = c_mean;
  };

  /**
   * Get the value c = beta/delta. This function may be useful when c
   * has been left to float, and the user wishes to fine out what is
   * can be reconstructed as.
   *
   * @return The value of c = beta/delta.
   */ 
  double get_c_mean(){
    return c_mean;
  };

  /** 
   * Set the constraint strength parameter for the magnitude.
   *
   * @param alpha1 The strength parameter for the magnitude.
   */
  void set_alpha1(double alpha1){
    this->alpha1 = alpha1;
  };

  /** 
   * Set the constraint strength parameter for the phase.
   *
   * @param alpha2 The strength parameter for the phase.
   */
  void set_alpha2(double alpha2){
    this->alpha2 = alpha2;
  };

  /** 
   * Get a new magnitude (amplitude) value calculated according to:
   *
   * \f[ A' = exp[ (1-\alpha_1)lnA + \alpha_1 \bar{C} \phi ] \f]
   *
   * Where A' is the new magnitude, A is the old magnitude, \f$ \phi
   * \f$ is the old phase, \f$\alpha_1\f$ controls the strength of the
   * constraint and \f$ \bar{C} \f$ is either the mean beta/delta, as
   * calculated in the parent TransmissionConstraint object, or the
   * fixed value passed by the use.
   *
   * @return The new magnitude of the transmission function.
   */
  double get_new_mag(double old_mag, double old_phase){
    return exp((1-alpha1)*log(old_mag) + alpha1*c_mean*old_phase);
  };

  /** 
   * Get a new phase value according to:
   *
   * \f[ 
   * \phi' = (1-\alpha_2)\phi + \alpha_2 \bar{C^{-1}} \ln A   
   * \f]
   *
   * Where \f$ \phi' \f$ is the new phase, A is the old magnitude, \f$
   * \phi \f$ is the old phase, \f$\alpha_2\f$ controls the strength of
   * the constraint and \f$ \bar{C^{-1}} \f$ is the inverse of either
   * the mean beta/delta, as calculated in the parent
   * TransmissionConstraint object, or the fixed value passed by the
   * use.
   *
   * @return The new phase of the transmission function.
   */
  double get_new_phase(double old_mag, double old_phase){
    return (1-alpha2)*old_phase + alpha2*log(old_mag)/c_mean;
  };

  /**
   *  Get the array which marks which elements of the transmission
   *  function array should have this constraint applied. Values of 0
   *  or less are excluded from the constraint. Positive values will
   *  have the constraint applied. This function is used by the
   *  TransmissionConstraint object to which it's given.
   *
   * @return A pointer to the array which maps which regions this
   * constraint should be applied to.
   */
  Double_2D * get_region(){
    return region;
  };

};



/**
 * @class TransmissionConstraint
 *
 * @brief An instance of this class will contain a description of any
 * constraints which should be applied to the transmission function
 * (FresnelCDI) or exit-surface-wave (PlanarCDI).
 *
 * This class allows constraints on the transmission function
 * (FresnelCDI) or exit-surface-wave (PlanarCDI), such as complex
 * constraints to be used. To achieve this the user needs to:
 * - create an instance of this class
 * - configure it (for example by setting some ComplexConstraint objects
 * - and finally by passing this object to the reconstruction (through 
 * the PlanarCDI or FresnelCDI function "set_complex_constraint").
 *
 * @see the example "ComplexConstraints_example.c".
 *
 */
class TransmissionConstraint{

 protected:

  /** a list of the complex_constraints */
  std::vector<ComplexConstraint*> complex_constraint_list;

  /** an array which hold the mapping between pixel and
      complex constraint */
  Double_2D * region_map;

  /** flag for whether unity on the transmission magnitude
      should be enforced */
  bool do_enforce_unity;

  /** flag for whether phase flipping should be enforced */
  bool do_charge_flip;

  /** the sign of the phase (+1, -1) to be flipped */
  int flip_sign;

  /** a function pointer to a cumstomized constraint */
  void (*custom_constraint)(Complex_2D&); 

 public:  
  
  /** The constructor, no parameters need to be passed */
  TransmissionConstraint();

  /** The destructor */
  ~TransmissionConstraint();
  

 /** 
  * delete all the complex constraints associated with
  * this TransmissionConstraint. Use with care! 
  **/
  void delete_complex_constraint_regions(){
    
    while(!complex_constraint_list.empty()) {
      delete complex_constraint_list.back();
      complex_constraint_list.pop_back();
    }
    
    delete region_map;
    region_map = 0;
    
  }


  /**
   * A complex constraint. See ComplexConstraint for details of this
   * class. Multiple regions can be defined by creating multiple
   * ComplexConstraint objects, and by adding them one at a time using
   * this function. All properties (other than the region) can be
   * reset at any point during the reconstruction (i.e. after calling
   * this function)
   *
   * @param new_constraint The complex constraint to add.
   */
  virtual void add_complex_constraint(ComplexConstraint & new_constraint);


  /**
   * Set whether charging flipping should or shouldn't be performed.
   * i.e. unphysical phases will be changed from phi -> -phi.  By
   * default it will be performed. This function also allows the user
   * to set the sign of the phases which should be flipped.
   *
   * @param enable true - flipping is performed, false - it is not.
   * @param flip_sign The sign of the phase to flip. If the value
   * passed is positive, all positive phases will be flipped to
   * negative. If a negative value is passed, the converse is true.
   */
  void set_charge_flipping(bool enable, int flip_sign=1){
    do_charge_flip = enable;
    if(flip_sign>0)
      this->flip_sign=1;
    else
      this->flip_sign=-1;
  };

  /**
   * Sets whether to enforce a unity constraint on the transmission
   * magnitude. If true, magnitudes of the transmission function which
   * are unphysical (above 1) will be reset to 1.
   * 
   * @param enable true - this constraint is applied, false - it is not.
   */
  void set_enforce_unity(bool enable){
    do_enforce_unity = enable;
  };
  
  /**
   * It's not possible to write code for every situation in which
   * someone wants to apply a constraint on the transmission function
   * or exit-surface-wave. For this reason, this function is provided.
   * It allows users to write their own custom constraint function and
   * have it executed without the need to recompiled the main
   * libraries of the software. The user's function must be of the
   * form: void function_name(Complex_2D & array). Where 'array' is
   * either the transmission function (FresnelCDI) or the
   * exit-surface-wave (PlanarCDI). 'array' is passed to the function,
   * and is also the return object (i.e. it is expected that the
   * function will alter 'array').
   * 
   * It will be executed after the default constraints. To only have
   * the user's function executed and no others, set this function and
   * set set_charge_flipping(false) and set_enforce_unity(false).
   * 
   * @param custom_constraint A function pointer (the name of the function)
   */
  void set_custom_constraint(void (*custom_constraint)(Complex_2D & transmission)){
    this->custom_constraint = custom_constraint;
  };

  
  /**
   * This is the function which is executed each time the support
   * constraint is applied in Fresnel or Planar CDI reconstruction.
   * The user should not generally need to use this function. An
   * exception would be if no constraint is applied during
   * reconstruction, but some general clean-up is required on the
   * final Transmission function. For example, applying the unity
   * constraint.
   *
   * Before applying any constraints, the mean values of c=beta/delta
   * will be calculated for each region defined through
   * ComplexConstraint objects. The mean is calculated as:
   *
   * \f[ 
   * \bar{C} = \frac{\sum{lnA}}{\sum{\phi}}
   * \f]
   *
   * Where A is the transmission function amplitude, \f$ \phi \f$ is
   * the phase and the summation runs over all values in the
   * transmission function array for the defined region.
   *
   * The constraints will be executed in the following order: - any
   * complex constraints which have been passed to the object - charge
   * flipping is performed - unity of the transmission function is
   * enforced - custom constraints are applied.
   */
  virtual void apply_constraint(Complex_2D & transmission);

};


#endif

