#include <iostream>
#include <string>
#include <cstdlib> 
#include <cmath>
#include "Complex_2D.h"
#include "PlanarCDI.h"
#include "io.h" 

using namespace std;

map<string,int> * PlanarCDI::algNameMap = PlanarCDI::set_up_algorithm_name_map();

map<string,int> * PlanarCDI::set_up_algorithm_name_map(){
  
  map<string,int> * temp_map = new map<string,int>;
  temp_map->insert(pair<string,int>("ER",ER));
  temp_map->insert(pair<string,int>("BIO",BIO));
  temp_map->insert(pair<string,int>("BOO",BOO));
  temp_map->insert(pair<string,int>("HIO",HIO));
  temp_map->insert(pair<string,int>("DM",DM));
  temp_map->insert(pair<string,int>("SF",SF));
  temp_map->insert(pair<string,int>("ASR",ASR));
  temp_map->insert(pair<string,int>("HPR",HPR));
  temp_map->insert(pair<string,int>("RAAR",RAAR));
  return temp_map;
}



/***************************************************************/
PlanarCDI::PlanarCDI(Complex_2D & initial_guess, int n_best)
  : complex(initial_guess),
    nx(initial_guess.get_size_x()),
    ny(initial_guess.get_size_y()),
    //    fft(nx,ny),
    temp_complex_PFS(nx,ny),
    temp_complex_PF(nx,ny),
    temp_complex_PS(nx,ny),
    temp_complex_PSF(nx,ny),
    beta(0.9),
    support(nx,ny),
    intensity_sqrt(nx,ny),
    n_best(n_best){
  
  set_algorithm(HIO);
  
  //initialize the best estimates
  best_array = new Complex_2D*[n_best];
  best_error_array = new double[n_best];
  for(int i=0; i<n_best; i++){
    best_array[i] = new Complex_2D(nx,ny);
    best_error_array[i]=1;
  }


  //initialize the beam-stop mask to null (not in use)
  beam_stop=0;

  //initialize the complex contraint to null;
  custom_complex_contraint = NULL; 

};

PlanarCDI::~PlanarCDI(){
  for(int i=0; i<n_best; i++)
    delete best_array[i];
  delete [] best_array;

}

Complex_2D * PlanarCDI::get_best_result(double & error, int index){
  if(index >=0 && index < n_best ){
    error = best_error_array[index];
    return best_array[index];
  }
  
  cout << "index out of bounds in PlanarCDI::get_best_result"<<endl;
  return NULL;
  
};

//double ** PlanarCDI::get_intensity_autocorrelation(){
void PlanarCDI::get_intensity_autocorrelation(Double_2D & autoc){

  //make a temporary object
  Complex_2D temp_intensity(nx,ny);

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      //set the real and imaginary components using the magnitude.
      double component = (1.0/sqrt(2.0))*(intensity_sqrt.get(i,j)*intensity_sqrt.get(i,j));
      temp_intensity.set_value(i,j,REAL, component);
      temp_intensity.set_value(i,j,IMAG, component);
    }
  }
  
  // fourier transform the intensity 
  //propagate_from_detector(temp_intensity);  

  temp_intensity.perform_backward_fft(); 
  temp_intensity.invert(true);

  //get the magnitude of the fourier transformed data.
  temp_intensity.get_2d(MAG, autoc);

  
}


void PlanarCDI::initialise_estimate(int seed){//change magnitude of guess HENRY
  //initialise the random number generator
  srand(seed);

  Complex_2D IGDP(nx,ny);

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
   
	IGDP.set_mag( i,j,intensity_sqrt.get(i,j) );

        double phase = 0.2*2*M_PI*(0.5-rand()/((double) RAND_MAX) );//initial value 0.2*pi
        IGDP.set_phase(i,j,phase);
     
      }
  }
  propagate_from_detector(IGDP);
  complex.copy(IGDP);
  
    for(int i=0; i<nx; i++){
      for(int j=0; j<ny; j++){

      if(!support.get(i,j)){ //enforce the support condition on the inital guess
	complex.set_value(i,j,REAL,0); 
	complex.set_value(i,j,IMAG,0);
      }
      else{
	   double r = support.get(i,j)*complex.get_value(i,j,REAL);
           double im = support.get(i,j)*complex.get_value(i,j,IMAG);
         
      	   complex.set_value(i,j,REAL,r);
           complex.set_value(i,j,IMAG,im);
      }
    }
  }




/*
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      if(!support.get(i,j)){ //enforce the support condition on the inital guess
	complex.set_value(i,j,REAL,0); 
	complex.set_value(i,j,IMAG,0);
      }
      else{
	double r = support.get(i,j)*( 255*(1/sqrt(2))*rand()/(double) RAND_MAX);//magnitude in r and im added by Henry
	double im = support.get(i,j)*( 255*(1/sqrt(2))*rand()/(double) RAND_MAX);
	
	complex.set_value(i,j,REAL,r); 
	complex.set_value(i,j,IMAG,im);
      }
    }
  }
*/
}

void PlanarCDI::set_support(const Double_2D & object_support, bool soften){
  double max = object_support.get_max();
  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){
      support.set(i,j,object_support.get(i,j)/max);
    }
  }
  //if required, convolve the support with a gaussian to soften the
  //edges
  if(soften) 
    convolve(support,3,5);

}

void PlanarCDI::set_beam_stop(const Double_2D & beam_stop_region){
  if(beam_stop==0)
    beam_stop = new Double_2D(nx, ny);
  beam_stop->copy(beam_stop_region);
}

void PlanarCDI::set_intensity(const Double_2D &detector_intensity){
  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){
      intensity_sqrt.set(i,j,sqrt(detector_intensity.get(i,j)));
    }
  }
}

double PlanarCDI::get_error(){
  return current_error;
}

void PlanarCDI::apply_support(Complex_2D & c){
  double support_value;

  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){

      support_value = support.get(i,j);
      if(support_value == 0){
	c.set_real(i,j,0);
	c.set_imag(i,j,0);
      }
      else if(support_value < 1){
	c.set_real(i,j,c.get_real(i,j) * support_value);
	c.set_imag(i,j,c.get_imag(i,j) * support_value);
      }
    }
  }

  if(custom_complex_contraint)
    (*custom_complex_contraint)(c);

}

void PlanarCDI::project_intensity(Complex_2D & c){
  propagate_to_detector(c);
  scale_intensity(c);
  propagate_from_detector(c);
}

void PlanarCDI::propagate_to_detector(Complex_2D & c){
  c.perform_forward_fft();
  c.invert(true);
}

void PlanarCDI::propagate_from_detector(Complex_2D & c){
  c.invert(true);
  c.perform_backward_fft(); 
}

void PlanarCDI::scale_intensity(Complex_2D & c){
  double norm2_mag=0;
  double norm2_diff=0;
  double current_int_sqrt=0;
  double current_mag=0;

  for(int i=0; i< nx; i++){
    for(int j=0; j< ny; j++){

      //reset the magnitude
      if(beam_stop==0 || beam_stop->get(i,j)>0){
	
	current_int_sqrt=intensity_sqrt.get(i,j);
	current_mag=c.get_mag(i,j);
	
	c.set_mag(i,j,current_int_sqrt);
      
	//calculate the error
	norm2_mag += current_int_sqrt*current_int_sqrt;
	norm2_diff += (current_mag-current_int_sqrt)
	  *(current_mag-current_int_sqrt);
      }
    }
  }
  current_error = (norm2_diff/norm2_mag);
  
}


void PlanarCDI::set_algorithm(int alg){
  /**  if(algorithm==alg){
    cout << "Warning you are trying to set the algorithm"
	 << " to the one already in use" << endl; 
    return;
    }**/

  algorithm = alg;
  
  switch(alg){

  case(ER): //checked
    set_custom_algorithm(0,0,0,1,0,0,0,0,0,0);
    break;
  case(BIO):
    set_custom_algorithm(0,beta,0,0,0,0,0,0,0,0);
    break;
  case(BOO):
    set_custom_algorithm(0,beta,0,0,0,0,0,0,1,0);
    break;
  case(HIO): //checked
    set_custom_algorithm(0,beta,1,0,0,0,0,0,0,0);
    break;
  case(DM): //checked
    set_custom_algorithm(beta,0,-1,0,-1,0,0,0,0,0);
    break;
  case(SF):
    set_custom_algorithm(0,1,0,1,0,0,0,0,0,0);
    break;
  case(ASR):
    set_custom_algorithm(0,1,1,0,0,0,0,0,0,0);
    break;
  case(HPR):
    set_custom_algorithm(0,beta,1,0,0,0,0,0,0,0);
    break;
  case(RAAR):
    set_custom_algorithm(0,1,beta,0,0,0,0,0,(1-beta),0);
    break;
  default:
    cout << "Algorithm unknown" <<endl;
  }

  //  print_algorithm();
}


void PlanarCDI::set_custom_algorithm(double m1, double m2, double m3, 
				      double m4, double m5, double m6, 
				      double m7, double m8,
				      double m9, double m10){
 
  algorithm_structure[PSF]=  m1 + m2 + m3 + m4;
  algorithm_structure[PFS]= -m1 + m5 + m6 + m7;
  algorithm_structure[PS] = -m3 - m6 - m8 + m10;
  algorithm_structure[PF] = -m2 - m5 + m8 + m9;
  algorithm_structure[PI] = -m4 - m7 - m9 - m10;
  
}

void PlanarCDI::print_algorithm(){
 
  cout << "Currently using the algorithm: "
       << "x(k+1) = x(k) + ("
       << algorithm_structure[PSF]<<"*PsPf + "
       << algorithm_structure[PFS]<<"*PfPs + "
       << algorithm_structure[PS]<<"*Ps + "
       << algorithm_structure[PF]<<"*Pf + "
       << algorithm_structure[PI]<<"*I"
       << ")x(k)"<< endl;
  cout << "Ps - support constraint, Pf - modulus constraint" << endl;
  
}

int PlanarCDI::iterate(){

  //PFS
  if(algorithm_structure[PFS]!=0){
    temp_complex_PFS.copy(complex);
    apply_support(temp_complex_PFS);
    project_intensity(temp_complex_PFS);
  }

  //F
  if(algorithm_structure[PF]!=0){
    temp_complex_PF.copy(complex);
    project_intensity(temp_complex_PF);
  }
  
  //S
  if(algorithm_structure[PS]!=0){
    temp_complex_PS.copy(complex);
    apply_support(temp_complex_PS);
  } 

  //SF
  if(algorithm_structure[PSF]!=0){
    if(algorithm_structure[PF]!=0){
      temp_complex_PSF.copy(temp_complex_PF);
      apply_support(temp_complex_PSF);
    }
    else{
      temp_complex_PSF.copy(complex);
      project_intensity(temp_complex_PSF);
      apply_support(temp_complex_PSF);
    }
  }

  //combine the result of the seperate operators
  double value_real, value_imag;
  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
 
      //Add the identity
      value_real = (1+algorithm_structure[PI])*complex.get_real(i,j);
      value_imag = (1+algorithm_structure[PI])*complex.get_imag(i,j);

      //Add the component from the PfPs operator
      if(algorithm_structure[PFS]!=0){
	value_real+=algorithm_structure[PFS]*temp_complex_PFS.get_real(i,j);
	value_imag+=algorithm_structure[PFS]*temp_complex_PFS.get_imag(i,j);
      }

      //Add the component from the Pf operator
      if(algorithm_structure[PF]!=0){
	value_real+=algorithm_structure[PF]*temp_complex_PF.get_real(i,j);
	value_imag+=algorithm_structure[PF]*temp_complex_PF.get_imag(i,j);
      }

      //Add the component from the Ps operator
      if(algorithm_structure[PS]!=0){
	value_real+=algorithm_structure[PS]*temp_complex_PS.get_real(i,j);
	value_imag+=algorithm_structure[PS]*temp_complex_PS.get_imag(i,j);
      }
      
      //Add the component from the PsPf operator
      if(algorithm_structure[PSF]!=0){
	value_real+=algorithm_structure[PSF]*temp_complex_PSF.get_real(i,j);
	value_imag+=algorithm_structure[PSF]*temp_complex_PSF.get_imag(i,j);
      }
      
      complex.set_real(i,j,value_real);
      complex.set_imag(i,j,value_imag);
      
    }
  }

  // check whether this estimate is as good as the current best
  // this is a bit dodgy since we are actually storing the 
  // estimate just after the best one.
  int place = 0;
  for( ; place < n_best && current_error > best_error_array[place]; place++); 
  
  //we found a new best estimate
  if(n_best>0 && place < n_best){

    Complex_2D * temp_pointer = best_array[n_best-1];
    temp_pointer->copy(complex);

    for(int i=n_best-1; i>place; i--){
      best_error_array[i] = best_error_array[i-1];
      best_array[i] = best_array[i-1];
    }
    
    best_error_array[place] = current_error;
    best_array[place] = temp_pointer;      
  }
  
  return SUCCESS;
}

void PlanarCDI::get_support(Double_2D & object_support){
  for(int i=0; i < nx; i++)
    for(int j=0; j < ny; j++)
      object_support.set(i,j,support.get(i,j));
}


void PlanarCDI::apply_shrinkwrap(double gauss_width, double threshold){
  
  Double_2D recon(nx,ny);
  complex.get_2d(MAG,recon);

  //convolve
  convolve(recon,gauss_width);
  
  //threshold
  apply_threshold(recon,threshold);

  set_support(recon);

}


void PlanarCDI::convolve(Double_2D & array, double gauss_width, 
			 int pixel_cut_off){
    //to speed up computation we only convolve 
  //up to 4 pixels away from the gaussian peak

  //make a temporary array to hold the smeared image
  Double_2D temp_array(nx,ny);
  
  //make a temporary array to hold the gaussian distribution.
  Double_2D gauss_dist(pixel_cut_off+1, pixel_cut_off+1);
  for(int i=0; i <= pixel_cut_off; i++){
    for(int j=0; j <= pixel_cut_off; j++){
      double denom = 2.0*gauss_width*gauss_width;
      gauss_dist.set(i,j,exp(-1*(i*i+j*j)/denom ) );
    }
  }      

  //now do the convolution
  //this is messy. First loop over the elements of
  //the array which was given as input
  double new_value;
  for(int i=0; i < nx; i++){
    for(int j=0; j < ny; j++){
      
      //now loop over the colvoluted array (the one we want to make).
      //Calculate the contribution to each element in it.
      
      new_value = 0;
      
      for(int i2=i-pixel_cut_off; i2 <= i+pixel_cut_off; i2++){
	for(int j2=j-pixel_cut_off; j2 <= j+pixel_cut_off; j2++){
	  if(i2<nx && i2>=0 && j2>=0 && j2<ny){
	    new_value += array.get(i2,j2)*gauss_dist.get(fabs(i-i2),fabs(j-j2));
	  }
	}
      }
      temp_array.set(i,j,new_value); 
    }
  }

  array.copy(temp_array);

}

/** threshold is a % of the maximum */
void PlanarCDI::apply_threshold(Double_2D & array, 
				double threshold){
  
  //find the maximum
  double max = array.get_max();
  
  //apply the threshold
  for(int i=0; i < nx; i++){
    for(int j=0; j < nx; j++){
      if( array.get(i,j) < (threshold*max) )
	array.set(i,j,0.0);
      else
	array.set(i,j,1.0);
    }
  }
}

void PlanarCDI::set_fftw_type(int type){
  temp_complex_PFS.set_fftw_type(type);
  temp_complex_PF.set_fftw_type(type);
  complex.set_fftw_type(type);
}

//**FUNCTION ADDED BY HENRY*********************//
double PlanarCDI::em(Complex_2D &a_o,bool fma, Complex_2D &b_o, bool fmb, Complex_2D &A,Complex_2D &B){
 
	Complex_2D a(nx,ny);
	Complex_2D b(nx,ny);
	a.copy(a_o);
	b.copy(b_o);

  
         
 	
	 if( fma ){
          
	  propagate_to_detector( a );
	   
	  for(int i = 0; i < nx; i++ ){
		for(int j = 0; j< ny; j++ ){
			a.set_mag(i,j, intensity_sqrt.get(i,j) );
		}
	  }

	  propagate_from_detector( a );
          
	}
	
 	 A.copy(a);
         
	 if( fmb ){
	   
	   propagate_to_detector( b ); 	   
         
	   for(int i = 0; i < nx; i++ ){
		for(int j = 0; j< ny; j++ ){
			b.set_mag(i,j, intensity_sqrt.get(i,j) );
		}
	  }

	   
	  propagate_from_detector( b ); 
         
	 }
	 
	 B.copy(b);
 	        

	double e_sq_sum = 0;
	double e_m = 0;
	double sum_int = 0;
	
     	for(int i = 0; i<nx; i++){
	   for(int j =0; j<ny; j++){
	
            e_sq_sum += (a.get_mag(i,j) - b.get_mag(i,j))*(a.get_mag(i,j) - b.get_mag(i,j));
   	
	    sum_int += intensity_sqrt.get(i,j)*intensity_sqrt.get(i,j);
	   }
	}
	
	e_m = e_sq_sum/( sum_int );
	
	return e_m;
     
	
}	
//**************************************************//	

		
	
	
          





