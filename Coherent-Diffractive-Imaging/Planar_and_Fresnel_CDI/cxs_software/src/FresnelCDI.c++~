#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include "Complex_2D.h"
#include "Double_2D.h"
#include "FresnelCDI.h"
#include "io.h" //
#include <sstream>

using namespace std;

FresnelCDI::FresnelCDI(Complex_2D & initial_guess,
		       Complex_2D & white_field,
		       double beam_wavelength,
		       double focal_detector_length,
		       double focal_sample_length,
		       double pixel_size,
		       double normalisation,
		       int n_best)
  :PlanarCDI(initial_guess,n_best),
   wavelength(beam_wavelength),
   pixel_length(pixel_size),
   norm(normalisation),
   illumination(nx,ny),
   B_s(nx,ny),
   B_d(ny,ny){

  illumination.copy(white_field);
  set_experimental_parameters(wavelength,focal_detector_length,
			      focal_sample_length,pixel_length);
  set_algorithm(ER);
}

void FresnelCDI::set_experimental_parameters(double beam_wavelength,
					     double focal_detector_length,
					     double focal_sample_length,
					     double pixel_size){
  
  wavelength = beam_wavelength;
  pixel_length = pixel_size;

  double x_mid = (nx-1)/2;
  double y_mid = (ny-1)/2;

  double zfd = focal_detector_length;
  double zfs = focal_sample_length;
  double zsd = focal_detector_length - focal_sample_length;

  double zsd_ = 1/(1/zsd - 1/focal_detector_length);

  double scaling_x = beam_wavelength*zsd/(pixel_length*nx);
  double scaling_y = beam_wavelength*zsd/(pixel_length*ny);

  double rho_2_d;
  double phi_B_d;
  double phi_B_s;

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      rho_2_d = pow(pixel_length*(x_mid-i),2) + 
	pow(pixel_length*(y_mid-j),2);

      phi_B_d = (M_PI*rho_2_d)/(beam_wavelength)*((1/focal_detector_length)-(1/zsd));
      phi_B_s = -phi_B_d;

      B_s.set_real(i,j,cos(phi_B_s));
      B_s.set_imag(i,j,sin(phi_B_s));

      B_d.set_real(i,j,cos(phi_B_d));
      B_d.set_imag(i,j,sin(phi_B_d));

    }
  }

}

void FresnelCDI::initialise_estimate(int seed){
  //initialise the random number generator
  srand(seed);
  
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      
      //do a simple initialisation
      
      //make the magnitude the difference of the intensity and
      //white-field
      double amp = intensity_sqrt.get(i,j)-norm*illumination.get_mag(i,j);
      
      //perturb the phase a bit about zero to allow random starts.
      double phase = 0.2*2*M_PI*(0.5-rand()/((double) RAND_MAX) );
      complex.set_mag(i,j,amp);
      complex.set_phase(i,j,phase);
      
    }
  }

  //take the result to the sample plane and apply the support
  propagate_from_detector(complex);
  apply_support(complex);
  
}


void FresnelCDI::scale_intensity(Complex_2D & c){

  c.add(illumination,norm); //add the white field

  PlanarCDI::scale_intensity(c);

  c.add(illumination,-norm);//subtract the white field

}

void FresnelCDI::propagate_from_detector(Complex_2D & c){
  c.multiply(B_d);
  c.perform_backward_fft();
  c.invert(true);
}

void FresnelCDI::propagate_to_detector(Complex_2D & c){
  c.invert(true); 
  c.perform_forward_fft();
  c.multiply(B_s);
}


void FresnelCDI::get_transmission_function(Complex_2D & result,
					   bool inforce_unity_mag){

  Complex_2D temp_illum(nx,ny);
  temp_illum.copy(illumination);

  //get the illuminating wavefield in the sample plane
  propagate_from_detector(temp_illum);

  //divide the estimate by the illuminating wavefield
  //and add unity.
  for(int i=0; i<nx; i++){
    for(int j=0; j<nx; j++){
      
      if(norm!=0&&temp_illum.get_mag(i,j)!=0){
      
	double new_mag = complex.get_mag(i,j)/(norm*temp_illum.get_mag(i,j));
	double new_phi = complex.get_value(i,j,PHASE) 
	  - temp_illum.get_value(i,j,PHASE);
	
	result.set_real(i,j,new_mag*cos(new_phi)+1);
	result.set_imag(i,j,new_mag*sin(new_phi));
	
	if(inforce_unity_mag && result.get_mag(i,j) > 1)
	  result.set_mag(i,j,1);      
      }
      else{
	result.set_real(i,j,0);
	result.set_imag(i,j,0);

      }
    }
  }
}



//**FUNCTION ADDED BY HENRY*********//
void FresnelCDI::get_the_transmission_function(Complex_2D & ESW, Complex_2D & result,
					   bool inforce_unity_mag){

  Complex_2D temp_illum(nx,ny);
  temp_illum.copy(illumination);

  //get the illuminating wavefield in the sample plane
  propagate_from_detector(temp_illum);

  //divide the estimate by the illuminating wavefield
  //and add unity.
  for(int i=0; i<nx; i++){
    for(int j=0; j<nx; j++){
      
      if(norm!=0&&temp_illum.get_mag(i,j)!=0){
      
	double new_mag = ESW.get_mag(i,j)/(norm*temp_illum.get_mag(i,j));
	double new_phi = ESW.get_value(i,j,PHASE) 
	  - temp_illum.get_value(i,j,PHASE);
	
	result.set_real(i,j,new_mag*cos(new_phi)+1);
	result.set_imag(i,j,new_mag*sin(new_phi));
	
	if(inforce_unity_mag && result.get_mag(i,j) > 1)
	  result.set_mag(i,j,1);      
      }
      else{
	result.set_real(i,j,0);
	result.set_imag(i,j,0);

      }
    }
  }
}

//////////////////////////////////////////////////////////
//** COMPLEX ESW EM*********************//
double FresnelCDI::C_f_em(Complex_2D &a_o,bool fma, Complex_2D &b_o, bool fmb){
 
	Complex_2D a(nx,ny);
	Complex_2D b(nx,ny);
	a.copy(a_o);
	b.copy(b_o);

  
         
 	
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
	
         
	 if( fmb ){
	   
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
	 
 	        
	double e_sq_sum = 0;
	double e_m = 0;
	double sum_int = 0;
	
     	for(int i = 0; i<nx; i++){
	   for(int j =0; j<ny; j++){
	
            e_sq_sum += pow(a.get_real(i,j) - b.get_real(i,j),2) + pow(a.get_imag(i,j) - b.get_imag(i,j), 2);
   	
	    sum_int += intensity_sqrt.get(i,j)*intensity_sqrt.get(i,j);
	   }
	}
	
	e_m = e_sq_sum/( sum_int );
	
	return e_m;
     
	
}	
//**************************************************//	


//**COMPLEX TF EM*********************//
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
     
	
}	
//**************************************************//	
///////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////
//** MAGNITUDE ESW EM*********************//
double FresnelCDI::M_f_em(Complex_2D &a_o,bool fma, Complex_2D &b_o, bool fmb){
 
	Complex_2D a(nx,ny);
	Complex_2D b(nx,ny);
	a.copy(a_o);
	b.copy(b_o);

  
         
 	
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
	
         
	 if( fmb ){
	   
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
	 
 	        
	double e_sq_sum = 0;
	double e_m = 0;
	double sum_int = 0;
	
     	for(int i = 0; i<nx; i++){
	   for(int j =0; j<ny; j++){
	
            e_sq_sum += pow(a.get_mag(i,j) - b.get_mag(i,j),2);
   	
	    sum_int += intensity_sqrt.get(i,j)*intensity_sqrt.get(i,j);
	   }
	}
	
	e_m = e_sq_sum/( sum_int );
	
	return e_m;
     
	
}	
//**************************************************//	


//**MAGNITUDE TF EM*********************//
double FresnelCDI::M_O_F_T_EM(Complex_2D &a_o,bool fma, Complex_2D &b_o, bool fmb, bool b_or_O){
 
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
	
            e_sq_sum += pow(A.get_mag(i,j) - B.get_mag(i,j), 2);
   	
	    sum_int += intensity_sqrt.get(i,j)*intensity_sqrt.get(i,j);
	   }
	}
	
	e_m = e_sq_sum/( sum_int );
	
	return e_m;
     
	
}	
//**************************************************//	
///////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//** PHASE ESW EM*********************//
double FresnelCDI::P_f_em(Complex_2D &a_o,bool fma, Complex_2D &b_o, bool fmb){
 
	Complex_2D a(nx,ny);
	Complex_2D b(nx,ny);
	a.copy(a_o);
	b.copy(b_o);

  
         
 	
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
	
         
	 if( fmb ){
	   
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
	 
 	        
	double e_sq_sum = 0;
	double e_m = 0;
	double sum_int = 0;
	
     	for(int i = 0; i<nx; i++){
	   for(int j =0; j<ny; j++){
	
            e_sq_sum += pow(a.get_value(i,j,PHASE) - b.get_value(i,j,PHASE),2);
   	
	    sum_int += intensity_sqrt.get(i,j)*intensity_sqrt.get(i,j);
	   }
	}
	
	e_m = e_sq_sum/( sum_int );
	
	return e_m;
     
	
}	
//**************************************************//	


//**PHASE TF EM*********************//
double FresnelCDI::P_O_F_T_EM(Complex_2D &a_o,bool fma, Complex_2D &b_o, bool fmb, bool b_or_O){
 
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
	
            e_sq_sum += pow(A.get_value(i,j,PHASE) - B.get_value(i,j,PHASE), 2);
   	
	    sum_int += intensity_sqrt.get(i,j)*intensity_sqrt.get(i,j);
	   }
	}
	
	e_m = e_sq_sum/( sum_int );
	
	return e_m;
     
	
}	
//**************************************************//	

void FresnelCDI::app_modulus(Complex_2D &cu){
 
 propagate_to_detector( cu );
 
	  cu.add(illumination,norm);
	  	  
	  for(int i = 0; i < nx; i++ ){
		for(int j = 0; j< ny; j++ ){
			cu.set_mag(i,j, intensity_sqrt.get(i,j) );
		}
	  }

	  cu.add(illumination,-norm);	
		
	  propagate_from_detector( cu );

}
///////////////////////////////////////////////////////////





