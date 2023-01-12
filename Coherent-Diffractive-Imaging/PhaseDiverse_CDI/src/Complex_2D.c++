// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

#include <iostream>  
#include <Complex_2D.h>
#include <stdlib.h>
#include <cstring>

using namespace std;

Complex_2D::Complex_2D(int x_size, int y_size){

  //set the array size
  nx = x_size;
  ny = y_size;

  //allocate memory for the array
#ifndef DOUBLE_PRECISION
  array = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*nx*ny);
#else
  array = (fftw_complex*) fftw_malloc(sizeof(fftwf_complex)*nx*ny);
#endif

  //initalise the fftw plans to null (not created yet. We will
  //create them when needed to avoid unnecessary time overhead).
  f_forward = 0;
  f_backward = 0;
  fftw_type = FFTW_MEASURE;
}

Complex_2D::Complex_2D(const Complex_2D& object){

  //set the array size
  nx = object.get_size_x();
  ny = object.get_size_y();


#ifndef DOUBLE_PRECISION
  array = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*nx*ny);//object.get_size_x()*object.get_size_y());
#else
  array = (fftw_complex*) fftw_malloc(sizeof(fftwf_complex)*object.get_size_x()*object.get_size_y());
#endif

  //initalise the fftw plans to null (not created yet. We will
  //create them when needed to avoid unnecessary time overhead).
  f_forward = 0;
  f_backward = 0;
  fftw_type = FFTW_MEASURE;

  for(int i=0; i < object.get_size_x(); i++){
    for(int j=0; j < object.get_size_y(); j++){

      set_real(i, j, object.get_real(i, j));
      set_imag(i, j, object.get_imag(i, j));
    }
  }
}

Complex_2D::~Complex_2D(){

#ifndef DOUBLE_PRECISION
  //free the memory of the array.
  fftwf_free(array);

  //free the memory of the fftw plans (but
  //only if it was actually allocated).
  if(f_forward)
    fftwf_destroy_plan(f_forward);
  if(f_backward)
    fftwf_destroy_plan(f_backward);

#else
  //free the memory of the array.
  fftw_free(array);

  //free the memory of the fftw plans (but
  //only if it was actually allocated).
  if(f_forward)
    fftw_destroy_plan(f_forward);
  if(f_backward)
    fftw_destroy_plan(f_backward);

#endif
}

//set the value at positions x,y. See Complex_2D.h for more info.
void Complex_2D::set_value(int x, int y, int component, double value){

  if(check_bounds(x,y)==FAILURE){
    cout << "can not set value out of array bounds" << endl;
    exit(1);
  }

  switch(component){

  case REAL :
    set_real(x,y,value);
    break;
  case IMAG :
    set_imag(x,y,value);
    break;
  case MAG :
    set_mag(x,y,value);
    break;
  case PHASE :
    set_phase(x,y,value);
    break;
  default:
    cout << "Value type in Complex_2D::set_value is unknown: " 
      << component << ". Must be REAL, IMAG, MAG or PHASE" << endl;
    exit(1);
  }
}

//get the value at positions x,y. See Complex_2D.h for more info.
double Complex_2D::get_value(int x, int y, int type) const {
  //by default we check that the value is within the bounds of the
  //array, but this can be turned off for optimisation.
  if(check_bounds(x,y)==FAILURE){
    cout << "can not get value out of array bounds" << endl;
    exit(1);
  }
  switch(type){
  case MAG:
    return get_mag(x,y);
  case REAL:
    return get_real(x,y);
  case IMAG:
    return get_imag(x,y);
  case PHASE: //goes between -pi and pi 
    return get_phase(x,y);
  case MAG_SQ:
    return pow(get_mag(x,y),2); //the square of the magnitude
  default:
    cout << "value type in Complex_2D::get_value is unknown" << endl;
    exit(1);
  }
}

//like get() but we do it for the entire array not just a single value.
void Complex_2D::get_2d(int type, Double_2D & result) const {

  for(int i=0; i < nx; i++)
    for(int j=0; j < ny; j++){
      result.set(i,j,get_value(i,j,type));
    }
}


//scale all the values in the array by the given factor.
void Complex_2D::scale(double scale_factor){

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){ 
      array[i*ny + j][REAL]*=scale_factor;
      array[i*ny + j][IMAG]*=scale_factor;

    }
  }

}

//add another Complex_2D to this one.
void Complex_2D::add(Complex_2D & c2, double scale){

  if(nx!=c2.get_size_x() || ny!=c2.get_size_y()){
    cout << "in Complex_2D::add, the dimensions of the "
      "input Complex_2D do not match the dimensions of "
      "this Complex_2D object" << endl;
    exit(1);
  }

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      array[i*ny + j][REAL]+=scale*c2.get_real(i,j);
      array[i*ny + j][IMAG]+=scale*c2.get_imag(i,j);
    }
  }
}

//multiply another Complex_2D with this one.
void Complex_2D::multiply(Complex_2D & c2, double scale){

  if(nx!=c2.get_size_x() || ny!=c2.get_size_y()){
    cout << "in Complex_2D::multiply, the dimensions of the "
      "input Complex_2D do not match the dimensions of "
      "this Complex_2D object" << endl;
    exit(1);
  }

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      // values are multiplied in the usual way
      // if c1 = a + ib and c2 = d + ie
      // then the new c1 is:
      // c1 = (a*d - b*e) + i(a*e + b*d)
      double new_real = c2.get_real(i,j)*this->get_real(i,j)
	- c2.get_imag(i,j)*this->get_imag(i,j);
      double new_imag = c2.get_real(i,j)*this->get_imag(i,j)
	+ c2.get_imag(i,j)*this->get_real(i,j);

      //and set the values
      set_real(i,j,scale*new_real);
      set_imag(i,j,scale*new_imag);
    }
  }
}

double Complex_2D::get_norm() const {

  double norm_squared=0;

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      norm_squared += pow(get_imag(i,j),2)+pow(get_real(i,j),2);
    }
  }
  return sqrt(norm_squared);
}

void Complex_2D::conjugate() {

  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      set_imag(i,j,-1*get_imag(i,j));
    }
  }

}

//make a new complex 2d that has idential values to this one
Complex_2D * Complex_2D::clone() const {

  Complex_2D * new_complex = new Complex_2D(nx,ny);
  for(int i=0; i < nx; ++i){
    for(int j=0; j < ny; ++j){
      new_complex->set_real(i,j, get_real(i,j));
      new_complex->set_imag(i,j, get_imag(i,j));
    }
  } 

  return new_complex;
}


//copy another array
void Complex_2D::copy(const Complex_2D & c){

  //check the bounds
  if(c.get_size_x()!=get_size_x()||
      c.get_size_y()!=get_size_y()){
    cout << "Trying to copy an array with different dimensions... "
      << "exiting"<<endl;
    exit(1);
  }

  //copy

#ifndef DOUBLE_PRECISION
  std::memcpy(array,c.array,sizeof(fftwf_complex)*nx*ny);
#else
  std::memcpy(array,c.array,sizeof(fftw_complex)*nx*ny);
#endif
}


//invert (and scale if we want to).
void Complex_2D::invert(bool scale){

  int middle_x = nx/2;
  int middle_y = ny/2;

  double scale_factor = 1;
  if(scale)
    scale_factor = 1.0/(sqrt(nx*ny));

  if(nx%2==1 || ny%2==1)
    cout << "WARNING: The array dimensions are odd "
      << "but we have assumed they are even when inverting an "
      << "array after FFT. This will probably cause you issues..."
      << endl;

  for(int i=0; i < nx; ++i){
    for(int j=0; j < middle_y; ++j){

      int j_new = j+middle_y; 
      int i_new = i+middle_x; 

      if(i >=  middle_x)
	i_new = i_new - 2*middle_x;

      double temp_rl = get_real(i_new,j_new);
      double temp_im = get_imag(i_new,j_new);

      set_real(i_new,j_new,get_real(i,j)*scale_factor);
      set_imag(i_new,j_new,get_imag(i,j)*scale_factor);

      set_real(i,j,temp_rl*scale_factor);
      set_imag(i,j,temp_im*scale_factor);
    }
  }

}


int Complex_2D::check_bounds(int x, int y) const{

  if(x < 0 || x >= nx || y < 0 || y >=ny )
    return FAILURE;

  return SUCCESS;
}

Complex_2D Complex_2D::get_padded(int x_add, int y_add){

  Complex_2D padded(nx+2*x_add, ny+2*y_add);

  for(int i=0; i<x_add; i++){
    for(int j=0; j<ny+2*y_add; j++){
      padded.set_real(i, j, 0);
      padded.set_imag(i, j, 0);
    }
  }

  for(int i=nx+x_add; i<nx+2*x_add; i++){
    for(int j=0; j<ny+2*y_add; j++){
      padded.set_real(i, j, 0);
      padded.set_imag(i, j, 0);
    }
  }

  for(int i=x_add; i<nx; i++){
    for(int j=0; j<y_add; j++){
      padded.set_real(i, j, 0);
      padded.set_imag(i, j, 0);
    }
  }

  for(int i=x_add; i<nx; i++){
    for(int j=ny+y_add; j<ny+2*y_add; j++){
      padded.set_real(i, j, 0);
      padded.set_imag(i, j, 0);
    }
  }

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      padded.set_real(i+x_add, j+y_add, get_real(i, j));
      padded.set_imag(i+x_add, j+y_add, get_imag(i, j));
      //std::cout<<i<<", "<<j<<", "<<x_add<<"\n";
    } 
  } 

  return padded;
}

Complex_2D Complex_2D::get_unpadded(int x_add, int y_add){

  Complex_2D unpadded(nx-2*x_add, ny-2*y_add);

  for(int i=0; i<nx-2*x_add; i++){
    for(int j=0; j<ny-2*y_add; j++){
      unpadded.set_real(i, j, get_real(i+x_add, j+y_add));
      unpadded.set_imag(i, j, get_imag(i+x_add, j+y_add));
      //      std::cout<<"Goodbye\n";

    }
  }
  return unpadded;
}




void Complex_2D::initialise_fft(){
  //creating the plan will erase the content of the array
  //so we need to be a bit tricky here.

#ifndef DOUBLE_PRECISION
  //make a new array 
  fftwf_complex * tmp_array;
  tmp_array = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*nx*ny);

  //make the plans
  f_backward = fftwf_plan_dft_2d(nx, ny, tmp_array, tmp_array, 
      FFTW_BACKWARD, fftw_type);
  f_forward = fftwf_plan_dft_2d(nx, ny,tmp_array,tmp_array, 
      FFTW_FORWARD, fftw_type);

  //now copy the array contents into the tmp_array,
  //free the old memory and update the pointer.
  std::memcpy(tmp_array,array,sizeof(fftwf_complex)*nx*ny);
  fftwf_free(array);

#else
  fftw_complex * tmp_array;
  tmp_array = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny);

  //make the plans
  f_backward = fftw_plan_dft_2d(nx, ny, tmp_array, tmp_array, 
      FFTW_BACKWARD, fftw_type);
  f_forward = fftw_plan_dft_2d(nx, ny,tmp_array,tmp_array, 
      FFTW_FORWARD, fftw_type);

  //now copy the array contents into the tmp_array,
  //free the old memory and update the pointer.
  std::memcpy(tmp_array,array,sizeof(fftw_complex)*nx*ny);
  fftw_free(array);

#endif

  array = tmp_array;

}

void Complex_2D::perform_forward_fft(){

  //make a new forward fft plan if we haven't made one already.
  if(f_forward==0 )
    initialise_fft();  

#ifndef DOUBLE_PRECISION
  fftwf_execute(f_forward);
#else
  fftw_execute(f_forward);
#endif
}


void Complex_2D::perform_backward_fft(){

  //make a new backward fft plan if we haven't made one already.
  if(f_backward==0)
    initialise_fft();
#ifndef DOUBLE_PRECISION
  fftwf_execute(f_backward);
#else
  fftw_execute(f_backward);
#endif
}


//this object is fourier transformed and the result placed in 'result'
void Complex_2D::perform_backward_fft_real(Double_2D & result){

#ifndef DOUBLE_PRECISION
  fftwf_plan fftw;
  fftw = fftwf_plan_dft_c2r_2d(nx,ny,array,result.array,FFTW_ESTIMATE); 
  fftwf_execute(fftw);
  fftwf_destroy_plan(fftw);
#else
  fftw_plan fftw;
  fftw = fftw_plan_dft_c2r_2d(nx,ny,array,result.array,FFTW_ESTIMATE); 
  fftw_execute(fftw);
  fftw_destroy_plan(fftw);  
#endif
}

//'result' is fourier transformed and the result placed in this object
void  Complex_2D::perform_forward_fft_real(Double_2D & input){
#ifndef DOUBLE_PRECISION
  fftwf_plan fftw;
  fftw = fftwf_plan_dft_r2c_2d(nx,ny,input.array,array,FFTW_ESTIMATE); 
  fftwf_execute(fftw);
  fftwf_destroy_plan(fftw);
#else
  fftw_plan fftw;
  fftw = fftw_plan_dft_r2c_2d(nx,ny,input.array,array,FFTW_ESTIMATE); 
  fftw_execute(fftw);
  fftw_destroy_plan(fftw);
#endif    
}



