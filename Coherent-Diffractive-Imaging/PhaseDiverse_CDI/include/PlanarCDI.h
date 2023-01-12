// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

/**
 * @file PlanarCDI.h
 * @class PlanarCDI
 * @author Nadia Davidson <nadiamd@unimelb.edu.au> 
 *
 * @brief The class which performs planar CDI reconstruction.
 *
 * The fundamental class used for performing CDI
 */

#ifndef PLANARCDI_H
#define PLANARCDI_H

//#include <map>
//#include <string>
#include "BaseCDI.h"

class PlanarCDI : public BaseCDI{

 public:

 PlanarCDI(Complex_2D & complex, unsigned int n_best=0)
   :BaseCDI(complex,n_best){};
  
  /**
   * Get the autocorrelation function of the intensity data.
   *
   * @param autoc The resulting autocorrelation.
   */ 
  virtual void get_intensity_autocorrelation(Double_2D & autoc);

  virtual void initialise_estimate(int seed=0);

  virtual void propagate_to_detector(Complex_2D & c);

  virtual void propagate_from_detector(Complex_2D & c);

};


#endif
