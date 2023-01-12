// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

/**
 * @file PhaseDiverseCDI.h
 * @class PhaseDiverseCDI
 * @author  Nadia Davidson <nadiamd@unimelb.edu.au> 
 *
 * @brief  The class which performs phase diverse and ptychographic 
 * reconstruction.
 *
 * This class can be used to perform phase diverse or ptychographic
 * reconstruction of either Fresnel or plane-wave CDI data. Any number
 * of frames (also called local/single frames or probes in this
 * documentation) may be added to the reconstruction. In order to
 * perform a reconstruction, you will need to create either a new
 * FresnelCDI or PlanarCDI object for each of these 'local' datasets.
 * The FresnelCDI or PlanarCDI objects must then be passed to a
 * PhaseDiverseCDI object. Because the each sub-iteration (ie. each
 * iteration of a local frame) is performed using the code in
 * FresnelCDI/PlanarCDI, all the functionality available in these
 * classes is also available here. For example complex constraint can
 * be set, shrink-wrap can be used, the estimate can be initialised
 * using the results from a previous reconstruction etc. An example
 * can be found in the /examples directory, demonstrating how to use
 * PhaseDiverseCDI.
 *
 * The displacement in position between 'local' frames may be either
 * transverse to the beam direction, or longitudinal, along the beam
 * direction. The longitudinal position (for FresnelCDI) is set when
 * the FresnelCDI objects are initially constructed. Transverse
 * positions are set when the FresnelCDI/PlanarCDI objects are added to
 * the PhaseDiverseCDI object. This class allows the transverse
 * positions to be automatically adjusted during reconstruction using
 * the "adjust_positions" function.
 *
 * The code allows for several options in the type of reconstruction
 * which is done:
 * <ul>
 * <li> The reconstruction can be performed in series or parallel.
 *   <ul>
 *   <li> In the case of series (see for example the paper....), a 'local'
 *   frame will undergo one or more iterations, the result will be
 *   updated to a 'global' estimate of the sample, and this estimate
 *   will form the starting point for the next frame. For each call to
 *   the method "iterate()" this process is repeated until each local
 *   frame is used once. The algorithm can be described by: <br> \f$
 *   T_{k+1} = (1-\beta w^n)T_k + \beta w^n T^n_k\f$ <br> where
 *   \f$T_{k+1}\f$ is the updated global function, \f$T^n\f$ is the
 *   updated local function for the nth local frame. \f$\beta\f$ is
 *   the relaxation parameter and the weight is \f$ w^n(\rho) =
 *   \alpha^n (\frac{|T^n(\rho)|}{max|T^n(\rho)|} )^\gamma \f$ for
 *   Fresnel CDI or \f$w^n = \alpha^n\f$ for Plane-wave CDI. The
 *   weight is zero outside of the support.
 *
 *   <li> For a parallel reconstruction (see .....), each frame with
 *   independently undergo one of more iteration, the result from all
 *   frames will be merged to form a new estimate of the sample, this
 *   estimate then becomes the starting point for the next iteration
 *   of all frames. The algorithm can be described by:
 *   <br> \f$ T_{k+1} = (1-\beta)T_k + \beta \sum_n(w^n T^n_k)\f$
 *   <br> where \f$T_{k+1}\f$, \f$T^n\f$, and \f$\beta\f$ were defined
 *   earlier. The weight, w, is similar to that used for series
 *   reconstruction, but the weight it normalised such that 
 *   \f$ \sum_n w^n = 1 \f$.  i.e. 
 *   \f$ w^n_{parallel}= w^n_{series} / \sum_n w^n_{series} \f$
 *  </ul> 
 *
 * <li> The number of local iterations to perform before updating the
 *   result to the 'global' function can be set.
 *
 * <li> The feedback parameter, beta, may be set. This quantity is used
 *   to set how much of the previous 'global' sample function will be
 *   left after the next 'global' iteration.
 *
 * <li> The amplification factor, gamma, and the probe scaling, alpha,
 *   may also be see. These parameters control the weighting of one
 *   frame (and pixels within a frame) with respect to each
 *   other.
 * </ul>
 */

#ifndef PHASED_H
#define PHASED_H

#include <BaseCDI.h>
#include <vector>

//forward declarations
class Complex_2D;
class PhaseDiverseCDI{

 protected:

  /* A vector holding a pointer to each of the FresnelCDI/PlanarCDI objects */
  std::vector<BaseCDI*> singleCDI;

  /* A vector holding a pointer to the transmission function/esw
     result for each FresnelCDI/PlanarCDI objects */
  std::vector<Complex_2D*> single_result;

  /** A vector holding the weighting function for each frame */  
  std::vector<Double_2D *> weights; 
  
  /* A vector of the transverse positions in x */
  std::vector<double> x_position;

  /* A vector of the transverse positions in y */
  std::vector<double> y_position;
  
  //parameters controlling the feedback
  double beta;
  double gamma;
  std::vector<double> alpha;

  /** The current estimate of the 'global' sample function */
  Complex_2D * object;

  /** The number of iterations to perform for each 'local' frame */
  int iterations_per_cycle;

  /** A factor which controls sub-pixel alignment. 
      1=regular, 2 = 2 'global' pixels for every 1 'local' pixel.
      i.e. 2 allows alignment to within half a pixel. */
  int scale;

  /** The size in x and y of the 'global' sample function. */
  int nx,ny;
  
  /** Minimum coordinates of the 'global' function based on the
   positions entered by the user. */
  int x_min;
  int y_min;

  /** a flag for running in either series or parallel mode */
  bool parallel; 

  /** a flag indicating whether the weights need to be recalculated */
  bool weights_set;

 public:

  enum {CROSS_CORRELATION,MINIMUM_ERROR};

  /** 
   * Construct a PhaseDiverseCDI object for phase diverse or
   * ptychographic reconstruction. The data should be entered later
   * using the 'add_new_position' function.
   *
   * @param beta The feedback parameter. By default this is 1 (no feedback).
   * @param gamma The amplification factor. By default this is 1 (no
   * amplification).
   * @param parallel true - run in parallel mode, false - run in series
   *        mode. By default series mode is set.
   * @param granularity  A factor which controls sub-pixel alignment. 
   *   1 = regular, 2 = 2 'global' pixels for every 1 'local' pixel.
   *   i.e. 2 allows alignment to within half a pixel. 
   *   NOTE: This is not currently working properly!
   */
  PhaseDiverseCDI(double beta=1.0, 
		  double gamma=1.0,
		  bool parallel=false,
		  int granularity=1);
  
  /** 
   * Destructor for PhaseDiverseCDI
   */
  ~PhaseDiverseCDI();


  /**
   *  
   *
   */
  void add_new_position(BaseCDI * local, 
			double x=0, double y=0, 
			double alpha=1);

  /**
   * Initialise the estimate of the 'global' object function. The
   * initialisation is performed using the current estimate for each
   * 'local' frame. For this reason you must first initialise each
   * FresnelCDI/PlanarCDI object individually before calling this
   * function.
   */
  void initialise_estimate();

  /** 
   * Perform one iteration. This with involve performing one or more
   * iterations for each sub-set of data. The exact implementation
   * depends on whether the reconstruction is being run in series or
   * parallel mode.
   */
  /*void*/double * iterate();
  
  /////////////////////////////////
  // Get and setter methods
  /////////////////////////////////

  /**
   * Set the number of 'local' frame iterations before updating the
   * result to the 'global' object function.
   *
   * @param iterations The number of 'local' iterations.
   */
  void set_iterations_per_cycle(int iterations){
    iterations_per_cycle = iterations;
  }

  /**
   * Set the feed-back parameter.
   *
   * @param beta The feedback parameter
   */
  void set_feedback_parameter(double beta){
    this->beta = beta;    
    weights_set = false;

  };

  void set_amplification_factor(double gamma){
    this->gamma = gamma;
    weights_set = false;
 
  };

  void set_probe_scaling(int n_probe, double alpha){
    if(this->alpha.size() > n_probe)
      this->alpha.at(n_probe)=alpha;
    else{
      std::cout << "In PhaseDiverseCDI::set_probe_scaling, "
	   << "the probe number given is too large. "
	   << "Are you trying to set the value of alpha "
		<< "before calling add_new_position?" <<std::endl;
    }
    weights_set = false;

  };
  
  /**
   * This function allows you to access the 'global' sample function.
   *
   * @return The current estimate of either the transmission (for
   * FresnelCDI) or exit-surface-wave (for PlanarCDI)
   */
  Complex_2D * get_transmission();
  
  /**
   * Set the 'global' sample function. This method could be used, for
   * example, instead of 'initialise_estimate' to initialise the
   * results to those from a previous reconstruction.
   *
   * @param new_transmission The transmission (for FresnelCDI) or
   * exit-surface-wave (for PlanarCDI) to copy.
   */
  void set_transmission(Complex_2D & new_transmission);
  
  ///////////////////////////////////////////
  // position adjustment
  //////////////////////////////////////////

  /**
   * Align the frames with respect to each other. The transverse positions will be 
   * automatically adjusted and the new positions used for the remainder of the 
   * reconstruction run. If needed, the new positions can be retrieved using the
   * "get_final_x(y)_position" functions. 
   *
   * Two methods to align the frames are provided:
   *
   *   - PhaseDiverse::CROSS_CORRELATION - implements the
   *     cross-correlation method using the magnitude of the
   *     transmission/exit-surface-wave.  There is pleanty of
   *     documentation describing how cross-correlation can be used
   *     for alignment. see http://en.wikipedia.org/wiki/Cross-correlation 
   *     for example. We use the fourier transform implementation.
   *     This method is useful if the relative positions are
   *     unknown to great than 10 pixels. It is also much faster than
   *     the alternative.
   *   
   *   - PhaseDiverse::MINIMUM_ERROR - aligns the frames by looking for the 
   *     position which gives the smallest error in the detector plane.
   *     The error metric is described in PlanarCDI.h. The algorithm works by
   *     checking the error in the 9 positions centered on the current positions: 
   *     (x-step_size, y-step_size), (x, y-step_size),  (x+step_size, y-step_size),
   *     (x-step_size, y) etc. Where "step_size" is an input parameter which
   *     has a default of 4 pixels. The position with the minimum error becomes 
   *     the new position, the step_size is decreased to half, and the procedure 
   *     is repeated. The algorithm ends when the step size is less than a pixel.
   *     This method of frame alignment is good for small refinement, but is less
   *     reliable if the relative differences are not known to within about 
   *     10-20 pixels. It is also significatly slower.
   * 
   * @param type either PhaseDiverse::CROSS_CORRELATION or 
   *             PhaseDiverse::MINIMUM_ERROR. PhaseDiverse::CROSS_CORRELATION
   *             is the default.
   * @param forwards By default the alignment is done in order of the frames, with the
   *                 second being aligned the the first, followed by the third to
   *                 the second and so on. In you needed to run the alignment in reverse 
   *                 e.g. aligning initially to the last frame, then this parameter 
   *                 should be false. You might use this, for example, if you want to check
   *                 that the alignment is consistent in both forward and backward directions.
   * @param x_min Positions below x+x_min pixels are not allowed. For 
   *                 PhaseDiverse::CROSS_CORRELATION, the algorithm will only search within 
   *                 the range x+x_min to x+x_max. For PhaseDiverse::MINIMUM_ERROR,
   *                 if a position below x_min is encountered, the algorithm aborts and
   *                 returns the positions to their original value.                
   * @param x_max See x_min.
   * @param y_min See y_min.
   * @param y_max See y_max.
   * @param step_size The initial step size to use in the case of PhaseDiverse::MINIMUM_ERROR.
   *
   */
  void adjust_positions(int type=CROSS_CORRELATION, 
			bool forwards=true,
			int x_min=-50, int x_max=+50,
			int y_min=-50, int y_max=+50,
			double step_size=4);

  /**
   * Get the current position of the 'n_probe'th frame.
   *
   * @param n_probe The index of the frame you wish to retrieve the
   *                position for. 0=the first frame you added, 
   *                1=2nd frame added etc....
   * @return the transverse horizontal position
   */
  double get_final_x_position(int n_probe){
    return x_position.at(n_probe);}
  
  /**
   * Get the current position of the 'n_probe'th frame.
   *
   * @param n_probe The index of the frame you wish to retrieve the
   *                position for. 0=the first frame you added, 
   *                1=2nd frame added etc....
   * @return the transverse vertical position
   */
  double get_final_y_position(int n_probe){
    return y_position.at(n_probe);}


 private:
  
  /**
   * Update a 'local' frame result to the 'global' object.
   * The result is weighted before being added.
   *
   * @param n_probe The local frame number.
   */
  void add_to_object(int n_probe);

  /**
   * Update a 'local' frame result from the 'global' object.
   *
   * @param n_probe The local frame number.
   */
  void update_from_object(int n_probe);

  /**
   * This function is basically a wrapper to check the BaseCDI type so
   * that the correct type of sample function is retrieve (either
   * transmission function or the exit-surface-wave.
   *
   * @param local A pointer to either a PlanarCDI or FresnelCDI object
   * @result result The estimate is copied from local to result
   */
  void get_result(BaseCDI * local, Complex_2D & result);

   /**
   * This function is basically a wrapper to check the BaseCDI type so
   * that the correct type of sample function is set (either
   * transmission function or the exit-surface-wave.
   *
   * @param local A pointer to either a PlanarCDI or FresnelCDI object
   * @result result The estimate is copied from result to local
   */ 
  void set_result(BaseCDI * local, Complex_2D & result);

  /**
   * This function is used to reallocate memory for the 'global'
   * sample object. It is only used when add_new_position is called,
   * in the case that the frame does not fit within the bounds of the
   * current object array.
   */
  void reallocate_object_memory(int new_nx, int new_ny);

  /**
   * This function is used during reconstruction in parallel mode. It
   * is similar to simply scaling all the elements of the 'global'
   * sample array, (e.g. with object.scale(factor). However, this
   * method preserves the value of the elements outside the boundary
   * of the sample. ie. it checks the weights of all the local frames
   * to work out which area to scale.
   * 
   * @param factor The scaled elements go from 'a' to 'factor*a'.
   *
   */
  void scale_object(double factor);
  
  /**  void get_object_sub_grid(Complex_2D & result,
			   double x_offset,
			   double y_offset); **/
 

  /**
   * Calculate the weights for each frame in the reconstruction.  See
   * a description above for how they are set.
   */
  void set_up_weights();

  /**
   * A position alignment function. See "adjust_positions" above for more detail.
   *
   * @param n_probe The index to the local frame which needs to be aligned
   * @param step_size Check the 8 positions located around the current positions using a
   *                  step size of "step_size" pixels.
   * @param min_x The minimum distance horizontally, in pixel, for which the re-alignment
   *              will check.
   * @param max_x See min_x.
   * @param min_y See min_y.
   * @param max_y See max_y.
   * @param tries How many times has this function been called for this positions. After 10
   *              tries the algorithm gives up and returns to the original position.
   */
  int check_position(int n_probe, double step_size=4, 
		     int min_x=-50, int max_x=50,
		     int min_y=-50, int max_y=50,
		     int tries = 0);
    
  /**
   * Get the global pixel "x" coordinate using a local frame "x" and
   * the frame offset.
   *
   * @param x the local frame pixel position in x 
   * @param x_offset the offset of the local frame 
   *                 with respect to the other frames.
   */
  inline int get_global_x_pos(int x, double x_offset){
    return x-x_offset-x_min; 
  }
 
  /**
   * See get_global_x_pos
   */
  inline int get_global_y_pos(int y, double y_offset){
    return y-y_offset-y_min;
  }

  /**
   * See get_global_x_pos
   */
  inline int get_local_x_pos(int x, double x_offset){
    return x + x_offset + x_min;
  }
  
  /**
   * See get_global_x_pos
   */
  inline int get_local_y_pos(int y, double y_offset){
    return y + y_offset + y_min; 
  }



};


#endif
