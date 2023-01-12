/**
 * @file Double_2D.h
 * @class Double_2D
 * @author  Nadia Davidson 
 * @date Last modified on 6/1/2011
 *
 * @brief A 2-dimensional array of doubles
 *
 * This class represents a 2D field of doubles. Setter and getter
 * methods are provided along with some other useful functions. This
 * is a convenient form to store and pass image data in the
 * reconstruction software.  A better implementation of this class
 * would use templates, but restricting the type to double is
 * sufficient our needs.
 */

#ifndef DOUBLE_2D_H
#define DOUBLE_2D_H

#include "string.h"

class Double_2D{

  /** the underlying 2-D array */
  double * array;
  
  /** the size in x */
  int nx;

  /** the size in y */
  int ny;

 public:

  /**
   * A constructor which creates an empty array (of no size).  Note
   * that memory has not been allocated if this method is used.
   */
  Double_2D();

  /**
   * Constructor that creates a 2D object with the given dimensions.
   * 
   * @param x_size The number of samplings in the horizontal direction
   * @param y_size The number of samplings in the vertical direction
   */
  Double_2D(int x_size, int y_size);
  
  /**
   * Destructor. Memory is deallocated here.
   */
  ~Double_2D();
     
  /**
   * Allocate memory for the array. This should only be used if
   * the constructor was called with no parameters!
   * 
   * @param x_size The number of samplings in the horizontal direction
   * @param y_size The number of samplings in the vertical direction
   */ 
  void allocate_memory(int x_size, int y_size);
  

  /**
   * Copy the contents of another Double_2D array to this one.  Note
   * that this is a quick copy method and no bounds checking is done.
   * 
   * @param double_array The array to copy from
   */
  void copy(const Double_2D & double_array);
  

  /**
   * Set the value at positions (x,y) WARNING: no bound checking is
   * done!
   *
   * @param x The horizontal position 
   * @param y The vertical position
   * @param value The value to set
   */
  inline void set(int x, int y, double value){
    array[x*ny+y]=value;
  };


  /**
   * Get the value at positions (x,y) WARNING: no bound checking is
   * done!
   *
   * @param x The horizontal position 
   * @param y The vertical position
   * @return The value at (x,y)
   */
  inline double get(int x, int y) const {
    return array[x*ny+y];
  };

  /**
   * Get the size in x;
   * 
   * @return The number of horizontal points.
   *  
   */
  inline int get_size_x() const {
    return nx;
  };

  /**
   * Get the size in y;
   * 
   * @return The number of vertical points.
   *  
   */
  inline int get_size_y() const {
    return ny;
  };


  /**
   * Get the sum of all values in the array. This is useful to
   * determine normalisation values.
   * 
   * @return The sum of all values in the array
   *  
   */
  double get_sum() const;


  /**
   * Get the maximum of all values in the array. 
   * 
   * @return The maximum value in the array
   *  
   */
  double get_max() const;

  /**
   * Get the minimum of all values in the array. 
   * 
   * @return The minimum value in the array 
   */
  double get_min() const;


  void add(const Double_2D & other_array, double norm=1.0){
    for(int i=0; i< nx; i++)
      for(int j=0; j< ny; j++)
	array[i*ny+j]+=norm*other_array.get(i,j);
  }

  void scale(double scale_by){
    for(int i=0; i< nx; i++)
      for(int j=0; j< ny; j++)
	array[i*ny+j]*=scale_by;
  }


};

#endif
