#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <cstdlib> 
#include "Complex_2D.h"
#include "Double_2D.h"
#include "FresnelCDI.h"
#include "TransmissionConstraint.h"
#include "io.h" //
#include <sstream>
#include "utils.h"

using namespace std;

#define FORWARD  -1
#define BACKWARD +1

FresnelCDI::FresnelCDI(Complex_2D & initial_guess,
		       Complex_2D & white_field,
		       double beam_wavelength,
		       double focal_detector_length,
		       double focal_sample_length,
		       double pixel_size,
		       double normalisation,
		       int n_best)
  :BaseCDI(initial_guess,n_best),
   illumination(nx,ny),
   norm(normalisation),
   coefficient(nx/2,ny/2)
{

  illumination.copy(white_field);
  illumination.scale(norm);
  
  illumination_at_sample=0;
  transmission = 0;
  set_experimental_parameters(beam_wavelength,focal_detector_length,
			      focal_sample_length,pixel_size);

  set_algorithm(ER);

}

FresnelCDI::~FresnelCDI(){

  if(illumination_at_sample)
    delete illumination_at_sample;
  
  if(transmission)
    delete transmission;

}


void FresnelCDI::set_experimental_parameters(double beam_wavelength,
					     double focal_detector_length,
					     double focal_sample_length,
					     double pixel_size){
  
  wavelength = beam_wavelength;
  pixel_length = pixel_size;
  this->focal_detector_length = focal_detector_length;
  this->focal_sample_length = focal_sample_length;
    
  double x_mid = (nx-1)/2.0;
  double y_mid = (ny-1)/2.0;

  double zfd = focal_detector_length;
  double zfs = focal_sample_length;
  double zsd = focal_detector_length - focal_sample_length;

  double factor = pixel_length*pixel_length*M_PI/(beam_wavelength)*((1/zfd) - (1/zsd));

  double phi;

  for(int i=0; i<nx/2.0; i++){
    for(int j=0; j<ny/2.0; j++){

      phi= factor*((x_mid-i)*(x_mid-i) + (y_mid-j)*(y_mid-j)); 

      coefficient.set_real(i,j,cos(phi));
      coefficient.set_imag(i,j,sin(phi));

    }
  }

  if(!illumination_at_sample)
    illumination_at_sample = new Complex_2D(nx,ny);
  
  illumination_at_sample->copy(illumination);
  propagate_from_detector(*illumination_at_sample);

}


void FresnelCDI::multiply_factors(Complex_2D & c, int direction){
  
  double x_mid = (nx-1)/2.0;
  double y_mid = (ny-1)/2.0;

  double old_real, old_imag;
  double coef_real, coef_imag;
  int i_, j_;

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      old_real = c.get_real(i,j);
      old_imag = c.get_imag(i,j);

      int i_ = x_mid - fabs(x_mid - i);
      int j_ = y_mid - fabs(y_mid - j);

      coef_real = coefficient.get_real(i_,j_);
      coef_imag = direction*coefficient.get_imag(i_,j_);

      c.set_real(i,j,old_real*coef_real - old_imag*coef_imag);
      c.set_imag(i,j,old_imag*coef_real + old_real*coef_imag);

    }
  }

  
}


void FresnelCDI::auto_set_norm(){
  Double_2D temp(nx,ny);
  illumination.get_2d(MAG,temp);
  double wf_norm = temp.get_sum();
  double int_norm = intensity_sqrt.get_sum();
  
  set_norm(int_norm/wf_norm);
  
}

void FresnelCDI::set_norm(double new_normalisation){

  double old_normalisation = norm;
  norm = new_normalisation;
  illumination.scale(new_normalisation/old_normalisation);
  
  if(illumination_at_sample)
    illumination_at_sample->scale(new_normalisation/old_normalisation);
  
}



void FresnelCDI::initialise_estimate(int seed){

  //initialise the random number generator
  srand(seed);
  Complex_2D random_trans(nx,ny);
  /**  Double_2D temp(nx,ny);
  illumination_at_sample->get_2d(REAL,temp);
  double max = temp.get_max();**/

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      double real = 1-0.2*rand()/((double) RAND_MAX);
      double imag = 0; //0.2*M_PI*rand()/((double) RAND_MAX) - 0.2*M_PI;

      random_trans.set_real(i,j,real);
      random_trans.set_imag(i,j,imag);

      //complex.set_real(i,j,max*rand()/((double) RAND_MAX));
      //complex.set_imag(i,j,0);//max*rand()/((double) RAND_MAX));
    }
  }
  
  set_transmission_function(random_trans,&complex);
  apply_support(complex);

  /**  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      complex.set_phase(i,j,0.0);
    }
    }**/

      //do a simple initialisation
      /**    
      //make the magnitude the difference of the intensity and
      //white-field
      double amp = intensity_sqrt.get(i,j)-illumination.get_mag(i,j);
      
      //perterb the phase a bit about zero to allow random starts.
      double phase = 0.2*2*M_PI*(0.5-rand()/((double) RAND_MAX) );
      complex.set_mag(i,j,amp);
      complex.set_phase(i,j,phase);
      
    }
  }

  //take the result to the sample plane and apply the support
  propagate_from_detector(complex);
  apply_support(complex);**/
  
}

void FresnelCDI::apply_support(Complex_2D & c){
  
  support_constraint(c);
  
  if(transmission_constraint){
    if(!transmission)
      transmission = new Complex_2D(nx,ny);

    get_transmission_function(*transmission,&c);
    set_transmission_function(*transmission,&c);
  }

}

void FresnelCDI::scale_intensity(Complex_2D & c){

  //  Double_2D result(nx,ny);
  
  //c.get_2d(MAG_SQ,result);
  //write_image("before_mag.tiff",result);
  //c.get_2d(PHASE,result);
  //write_image("before_p.tiff",result);
  
  c.add(illumination); //add the white field

  //c.get_2d(REAL,result);
  //write_image("after_r.tiff",result);
  //c.get_2d(MAG_SQ,result);
  //write_image("after_mag.tiff",result);
  BaseCDI::scale_intensity(c);

  c.add(illumination,-1.0);//subtract the white field

}

void FresnelCDI::propagate_from_detector(Complex_2D & c){
  //  c.multiply(B_d);
  multiply_factors(c,BACKWARD);
  c.perform_backward_fft();
  c.invert(true);
}

void FresnelCDI::propagate_to_detector(Complex_2D & c){
  c.invert(true); 
  c.perform_forward_fft();
  //c.multiply(B_s);
  multiply_factors(c,FORWARD);
}


void FresnelCDI::set_transmission_function(Complex_2D & transmission,
					   Complex_2D * esw){
  if(!esw)
    esw=&complex;

  check_illumination_at_sample();
  /**  if(!illumination_at_sample){
    illumination_at_sample = new Complex_2D(nx,ny);
    illumination_at_sample->copy(illumination);
    propagate_from_detector(*illumination_at_sample);
    }**/

  double ill_r;
  double trans_r;
  double ill_i;
  double trans_i;
  
  for(int i=0; i<nx; i++){
    for(int j=0; j<nx; j++){
      ill_r = illumination_at_sample->get_real(i,j);
      trans_r = transmission.get_real(i,j);
      ill_i = illumination_at_sample->get_imag(i,j);
      trans_i = transmission.get_imag(i,j);
   
      // ESW = TL - L
      // T - transmission function, L - illumination
      esw->set_real(i,j,ill_r*trans_r - ill_i*trans_i - ill_r);
      esw->set_imag(i,j,ill_r*trans_i + ill_i*trans_r - ill_i);      
    }
  }
}

void FresnelCDI::check_illumination_at_sample(){ 
  if(!illumination_at_sample){
    illumination_at_sample = new Complex_2D(nx,ny);
    illumination_at_sample->copy(illumination);
    propagate_from_detector(*illumination_at_sample);
  }
}

const Complex_2D & FresnelCDI::get_illumination_at_sample(){ 
  check_illumination_at_sample();
  return *illumination_at_sample;
}

void FresnelCDI::get_transmission_function(Complex_2D & result, 
					   Complex_2D * esw){

  //divide the estimate by the illuminating wavefield and add unity.
 
  //the code below could written more eligantly with get_mag and get_phase etc.
  //but the code below is faster because it doesn't use the maths tan function.
  //This is important when complex constraints are applied

  if(!esw)
    esw=&complex;
  check_illumination_at_sample();

  double ill_r;
  double esw_r;
  double ill_i;
  double esw_i;
  
  //  double combined_r;
  //double combined_i;

  double real_numerator;
  double imag_numerator;
  double denom;
  
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
   
      ill_r = illumination_at_sample->get_real(i,j);
      ill_i = illumination_at_sample->get_imag(i,j);
      denom = ill_r*ill_r + ill_i*ill_i;
      
      if(denom!=0){
  
	esw_r = esw->get_real(i,j);
	esw_i = esw->get_imag(i,j);

	real_numerator = ill_r*(esw_r+ill_r) + ill_i*(esw_i+ill_i);
	imag_numerator = ill_r*(esw_i+ill_i) - ill_i*(esw_r+ill_r);

	result.set_real(i,j,real_numerator/denom);
	result.set_imag(i,j,imag_numerator/denom);
	
	
	//	if(inforce_unity_mag && result.get_mag(i,j) > 1)
	//	  result.set_mag(i,j,1);      
      }
      else{
	result.set_real(i,j,1.0);
	result.set_imag(i,j,0.0);

      }
    }
  }
  
  if(transmission_constraint){
    transmission_constraint->apply_constraint(result);
  }
}


//**COMPLEX TF EM*********************// 
/*
double FresnelCDI::C_O_F_T_EM(Complex_2D &a_o,bool fma, Complex_2D &b_o, bool fmb, bool b_or_O){
 
	Complex_2D a(nx,ny);
        Complex_2D b(nx,ny);
	a.copy(a_o);
        b.copy(b_o);
        Complex_2D A(nx,ny);
        Complex_2D B(nx,ny);

  
         
 	
	if( fma ){
          
	  propagate_to_detector( a );
 
	  a.add(illumination,norm);
	  	  
	  for(int i = 0; i < nx; i++ ){
		for(int j = 0; j< ny; j++ ){
			a.set_mag(i,j, intensity_sqrt.get(i,j) );
		}
	  }

	  a.add(illumination,-norm);	
		
	  propagate_from_detector( a );
          
	}
	
        get_the_transmission_function(a, A, false);

         
     if(b_or_O=1)
      {           
	if( fmb )
       {
          
	  propagate_to_detector( b );
 
	  b.add(illumination,norm);
	  	  
	  for(int i = 0; i < nx; i++ ){
		for(int j = 0; j< ny; j++ ){
			b.set_mag(i,j, intensity_sqrt.get(i,j) );
		}
	  }

	  b.add(illumination,-norm);	
		
	  propagate_from_detector( b );
          
	}
	
        get_the_transmission_function(b, B, false);
        
         
      } 
    
    else if(b_or_O=0){
      B.copy(b);} 

         
	
	double e_sq_sum = 0;
	double e_m = 0;
	double sum_int = 0;
	
     	for(int i = 0; i<nx; i++){
	   for(int j =0; j<ny; j++){
	
            e_sq_sum += pow(A.get_real(i,j) - B.get_real(i,j), 2)+ pow(A.get_imag(i,j) - B.get_imag(i,j), 2);
   	
	    sum_int += intensity_sqrt.get(i,j)*intensity_sqrt.get(i,j);
	   }
	}
	
	e_m = e_sq_sum/( sum_int );
	
	return e_m;
     
	
}*/	
//**************************************************//	
///////////////////////////////////////////////////////////



/**double FresnelCDI::refine_sample_to_focal_length(double min, double max,
				     int points_per_scan,
				     double precision,
				     int iterations_before_comparison,
				     int crop_min_x, int crop_min_y,
				     int crop_max_x, int crop_max_y ){

  //set some defaults
  if(min==0&&max==0){
    min = 0.8*focal_sample_length;
    max = 1.2*focal_sample_length;
  }
  if(precision==0)
    precision = 0.01*focal_sample_length;
  if(crop_max_x==0)
    crop_max_x=nx;
  if(crop_max_y==0)
    crop_max_y=ny;

  
  //make copies of the object esw estimate ("complex")
  //and support. The esw and support will be overriden
  //during this function call, and we will use the copied
  //to reset them back to their original state at the end.
  Complex_2D * complex_copy = new Complex_2D(nx,ny);
  complex_copy->copy(complex);
  Double_2D * support_copy = new Double_2D(nx,ny);
  support_copy->copy(get_support());
  //double fd = focal_detector_length;
  double fs = focal_sample_length;

  //store the focal metric values in this vector object.
  vector<double> focal_metrics;
  
  //test the images at several points between the min and max.
  //focal lengths.
  for(double f = min; f <= max; f+=(max-min)/points_per_scan){
    
    //work out the new length and set it.
    cout << endl << f << endl;
      
    //maintain the detector sample distance
    //    double fd_new = fd + fs - f;  //focal to detector
    double fs_new = f;//    ; //focal to sample

    //reset the initial estimate so that all lengths
    //are compared from the same starting point
    
    set_experimental_parameters(wavelength,
				focal_detector_length,
				fs,
				pixel_length);
    complex.copy(*complex_copy);
    propagate_to_detector(complex);
    set_experimental_parameters(wavelength,
				focal_detector_length,
				fs_new,
				pixel_length);
    propagate_from_detector(complex);
    
    //scale the support for size
    double pixel_ratio = (focal_detector_length-fs)/(focal_detector_length-fs_new);
    cout << "pixel_ratio="<<pixel_ratio<<endl; 
    double length_ratio = fs_new/fs;
    cout << "length_ratio="<<length_ratio<<endl;
    double scale = length_ratio/pixel_ratio;
    cout << "scale="<<scale<<endl;
    Double_2D temp_support(nx,ny);
    temp_support.copy(*support_copy);
    rescale(temp_support,1.0/scale);
    set_support(temp_support);
     
    //initialise_estimate();
    
    for(int i=0; i<iterations_before_comparison; i++) 
      iterate();

   
    Double_2D result(nx,ny);
    complex.get_2d(MAG_SQ,result);

    rescale(result,scale);
    result.scale(scale*scale);
    
    Double_2D crop_result(crop_max_x-crop_min_x,crop_max_x-crop_min_y);
    crop(result,crop_result,crop_min_x,crop_min_y);

    //temp remove later
    char buff[80];
    //     static int counter=0;
    sprintf(buff,"temp_%f.tiff",f);
    write_image(buff,result);
    //   counter++;
        
     //calculate the focal metric
     double metric_value = edges(crop_result);
     focal_metrics.push_back(metric_value);

  }

  double largest = focal_metrics.at(0);
  double smallest = focal_metrics.at(0);
  double sum = 0;
  int best_index = 0;

  //fine the largest value of the focal metric.

  for(int l=0; l<focal_metrics.size(); l++){
    if(largest<focal_metrics.at(l)){
      largest=focal_metrics.at(l);
      best_index = l;
    }
    if(smallest>focal_metrics.at(l))
      smallest=focal_metrics.at(l);
    sum+=focal_metrics.at(l);
    
    cout //<< "("<< (min + (max-min)*l/points_per_scan) 
      << "," << focal_metrics.at(l); // << ") ";
    
  }

  //return everything to normal
  set_support(*support_copy);
  complex.copy(*complex_copy);
  focal_sample_length = fs;
  delete support_copy;
  delete complex_copy;


  //if we have already acheive the required 
  //precision return as we're finished.
  if(max-min < precision){
    double new_focal_sample_length = min + (max-min)*best_index/points_per_scan;
    set_experimental_parameters(wavelength,
				focal_detector_length,
				new_focal_sample_length,
				pixel_length);

    return new_focal_sample_length;
  }

  double new_min = min;
  double new_max = max;
  //  double cut_off = 0.5*(largest-smallest)+smallest;
  double largest_below = min;
  double largest_above = max;

  //work out new min and max values
  for(int l=0; l<best_index; l++){
    if( focal_metrics.at(l) > largest_below ){
      largest_below = focal_metrics.at(l);
      new_min = min + (max-min)*l/points_per_scan;
    }
  }

  for(int l=focal_metrics.size()-1; l>best_index; l--){
    if( focal_metrics.at(l) > largest_above ){
      largest_above = focal_metrics.at(l);
      new_max = min + (max-min)*l/points_per_scan;
    }
  }
  
  if(new_min==min&&new_max==max){
    cout << "Could not refine the focal length. "
	 << "A maximum of the focal metric was not found "
	 << "because the function contained multiple local maxima. "
	 << "Restoring the original focal-sample length." 
	 << endl;
    return fs;
  }

  //reset the focal-sample length to the best found so far
  //focal_sample_length = min + (max-min)*best_index/points_per_scan;

  //recursively call this function until the error is small enough.

  //  int min_index = best_index - points_per_scan/4;
  //int max_index = best_index + points_per_scan/4;

  double best_length = refine_sample_to_focal_length(new_min, new_max,
						     points_per_scan,
						     precision,
						     iterations_before_comparison,
						     crop_min_x, crop_min_y,
						     crop_max_x, crop_max_y);
  


  //return the value of focal length found.
  return best_length;

  }**/
