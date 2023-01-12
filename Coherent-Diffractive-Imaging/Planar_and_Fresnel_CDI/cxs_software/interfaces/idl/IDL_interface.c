//standard headers
#include <cstdlib>
#include <iostream>
#include <typeinfo>
#include <sstream>
#include <string>

//idl header
#include "idl_export.h"

//cxs software headers
#include "Complex_2D.h"
#include "Double_2D.h"
#include "PlanarCDI.h"
#include "FresnelCDI.h"
#include "FresnelCDI_WF.h"
#include "io.h"

Complex_2D * esw = 0;
PlanarCDI * reco = 0;
int total_iters = 0;
using namespace std;

/************ helper methods ******************/

//methods for converting between the IDL types and the C++ types.
void copy_to_double_2d(Double_2D & cxs_array, double * IDL_array){

  int nx = cxs_array.get_size_x();
  int ny = cxs_array.get_size_y();

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      cxs_array.set(j,nx-1-i,*(IDL_array));
      IDL_array++;
    }
  }
}

void copy_from_double_2d(Double_2D & cxs_array, double * IDL_array){

  int nx = cxs_array.get_size_x();
  int ny = cxs_array.get_size_y();

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      *(IDL_array)= cxs_array.get(j,nx-1-i);
      IDL_array++;
    }
  }
}

void copy_to_complex_2d(Complex_2D & cxs_array, IDL_COMPLEX * IDL_array){

  int nx = cxs_array.get_size_x();
  int ny = cxs_array.get_size_y();

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      cxs_array.set_real(j,nx-1-i,(*IDL_array).r);
      cxs_array.set_imag(j,nx-1-i,(*IDL_array).i);
      IDL_array++;
    }
  }
}

void copy_from_complex_2d(Complex_2D & cxs_array, IDL_COMPLEX * IDL_array){

  int nx = cxs_array.get_size_x();
  int ny = cxs_array.get_size_y();  

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      (*IDL_array).r= cxs_array.get_real(j,nx-1-i);
      (*IDL_array).i= cxs_array.get_imag(j,nx-1-i);
      IDL_array++;
    }
  }
}

//error checking method. Make sure the objects have memory allocated.
void check_objects(){
  if(!esw || !reco){
    char buf[] = "You need to call one of the CXS_INIT_ routines before calling this routine";
    IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, buf); 
    exit(1);
  }
}


//Allows the dimensions of the image to be passed back to the IDL code
extern "C" IDL_LONG IDL_get_array_x_size(int argc, void *argv[]){
  if(!esw)
    return 0;
  else
    return esw->get_size_x();
}

extern "C" IDL_LONG IDL_get_array_y_size(int argc, void *argv[]){
  if(!esw)
    return 0;
  else
    return esw->get_size_y();
}





/************ input-output methods ******************/


extern "C" void IDL_read_dbin(int argc, void * argv[])
{
  int ny = *(int*) argv[0];
  int nx = *(int*) argv[1];

  IDL_STRING filename = *(IDL_STRING*)argv[2];
  Double_2D temp(nx,ny);
  read_dbin(filename.s,nx,ny,temp);
  copy_from_double_2d(temp, (double*) argv[3]);

}

extern "C" void IDL_read_ppm(int argc, void * argv[])
{
  int ny = *(int*) argv[0];
  int nx = *(int*) argv[1];

  IDL_STRING filename = *(IDL_STRING*)argv[2];
  Double_2D temp(nx,ny);
  read_ppm(filename.s,temp);
  copy_from_double_2d(temp, (double*) argv[3]);

}



extern "C" void IDL_read_cplx(int argc, void * argv[])
{
  int ny = *(int*) argv[0];
  int nx = *(int*) argv[1];

  IDL_STRING filename = *(IDL_STRING*)argv[2];
  Complex_2D temp(nx,ny);
  read_cplx(filename.s,temp);
  copy_from_complex_2d(temp, (IDL_COMPLEX*) argv[3]);

}

extern "C" void IDL_read_tiff(int argc, void * argv[])
{

  IDL_STRING filename = *(IDL_STRING*)argv[2];
  Double_2D temp;
  read_tiff(filename.s,temp);
  copy_from_double_2d(temp, (double*) argv[3]);

}

extern "C" void IDL_write_dbin(int argc, void * argv[])
{
  int ny = *(int*) argv[0];
  int nx = *(int*) argv[1];
  
  IDL_STRING filename = *(IDL_STRING*)argv[3];
  Double_2D temp(nx,ny);
  copy_to_double_2d(temp, (double*) argv[2]);
  write_dbin(filename.s,temp);
}

extern "C" void IDL_write_cplx(int argc, void * argv[])
{
  int ny = *(int*) argv[0];
  int nx = *(int*) argv[1];
  
  IDL_STRING filename = *(IDL_STRING*)argv[3];

  Complex_2D temp(nx,ny);

  copy_to_complex_2d(temp, (IDL_COMPLEX*) argv[2]);

  write_cplx(filename.s,temp);
}


/************ reconstuction methods ******************/



/******* memory allocation and deallocation methods ****/

extern "C" void IDL_deallocate_memory(int argc, void * argv[])
{
  if(reco!=0){
    delete reco ;
    reco = 0 ;
  }

  if(esw!=0){
    delete esw ;
    esw = 0 ;
  }
}


void common_init(int argc, void * argv[], int max_args){


  IDL_deallocate_memory(0,0);

  total_iters = 0;

  int ny = *(int*) argv[0];
  int nx = *(int*) argv[1];

  esw = new Complex_2D(nx,ny);  
  
  if(argc == max_args)
    copy_to_complex_2d(*esw,(IDL_COMPLEX*) argv[max_args-1]);

}


extern "C" void IDL_planar_init(int argc, void * argv[])
{
  common_init(argc, argv, 3);
  reco = new PlanarCDI(*esw);
}

extern "C" void IDL_fresnel_wf_init(int argc, void * argv[])
{
 
  common_init(argc, argv, 7); 
  reco = new FresnelCDI_WF(*esw, 
			   *(double*) argv[2],
			   *(double*) argv[3], 
			   *(double*) argv[4],
			   *(double*) argv[5]);
}


extern "C" void IDL_fresnel_init(int argc, void * argv[])
{

  common_init(argc, argv, 9); 

  Complex_2D white_field(*(int*)argv[1],*(int*)argv[0]);

  copy_to_complex_2d(white_field,(IDL_COMPLEX*) argv[2]);

  reco = new FresnelCDI(*esw,
			white_field,
			*(double*) argv[3],
			*(double*) argv[4], 
			*(double*) argv[5],
			*(double*) argv[6],
			*(double*) argv[7]);

}

/***** Getter and setter methods ********************/


extern "C" void IDL_set_support(int argc, void * argv[]){

  check_objects();

  int ny = *(int*) argv[0];
  int nx = *(int*) argv[1];

  Double_2D temp_2D(nx,ny);
  copy_to_double_2d(temp_2D,(double*) argv[2]); 

  reco->set_support(temp_2D);
}

extern "C" void IDL_set_beam_stop(int argc, void * argv[]){

  check_objects();

  int ny = *(int*) argv[0];
  int nx = *(int*) argv[1];

  Double_2D temp_2D(nx,ny);
  copy_to_double_2d(temp_2D,(double*) argv[2]); 

  reco->set_beam_stop(temp_2D);
}

extern "C" void IDL_set_intensity(int argc, void * argv[]){

  check_objects();

  int ny = *(int*) argv[0];
  int nx = *(int*) argv[1];

  Double_2D temp_2D(nx,ny);
  copy_to_double_2d(temp_2D,(double*) argv[2]); 

  reco->set_intensity(temp_2D);
}

extern "C" void IDL_initialise_esw(int argc, void * argv[]){
  check_objects();
  reco->initialise_estimate(*(IDL_LONG*) argv[0]);
}

extern "C" void IDL_iterate(int argc, void * argv[]){
  check_objects();
  int iterations = *(IDL_LONG*) argv[0];

  for(int i=0; i<iterations; i++){
    total_iters++;
    reco->iterate();

    ostringstream oss (ostringstream::out);
    oss << "Error for iteration " << total_iters - 1 
	<< " is "<< reco->get_error() << endl
	<< "Iteration: " << total_iters << endl;
    IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());    

  }

  copy_from_complex_2d(*esw,(IDL_COMPLEX*) argv[1]); 

}

extern "C" void IDL_set_algorithm(int argc, void * argv[]){
  check_objects();
  IDL_STRING alg_name = *(IDL_STRING*)argv[0];
  int alg = PlanarCDI::getAlgFromName(alg_name.s);
  reco->set_algorithm(alg);
}

extern "C" void IDL_set_custom_algorithm(int argc, void * argv[]){
  check_objects();
  reco->set_custom_algorithm(*(double*)argv[0],
			     *(double*)argv[1],
			     *(double*)argv[2],
			     *(double*)argv[3],
			     *(double*)argv[4],
			     *(double*)argv[5],
			     *(double*)argv[6],
			     *(double*)argv[7],
			     *(double*)argv[8],
			     *(double*)argv[9]);
  reco->print_algorithm();

}


extern "C" void IDL_set_relaxation_parameter(int argc, void * argv[]){
  check_objects();
  reco->set_relaxation_parameter(*(double*) argv[0]);
}

extern "C" void IDL_apply_shrinkwrap(int argc, void * argv[]){
  check_objects();
  double gauss_width = *(double*) argv[0];
  double threshold = *(double*) argv[1];

  
  ostringstream oss (ostringstream::out);
  oss << "Applying shrink wrap with a gaussian width of "
      << gauss_width <<" pixels and a threshold of " <<  threshold
      << "  the maximum pixel value" <<endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str()); 
  reco->apply_shrinkwrap(gauss_width,threshold);
}


extern "C" void IDL_get_best_result(int argc, void * argv[]){
  check_objects();
  Complex_2D * temp;
  double error;
  temp = reco->get_best_result(error);
  copy_from_complex_2d(*temp,(IDL_COMPLEX*) argv[0]);

  ostringstream oss (ostringstream::out);
  oss << "Error for the best result so far is "<< error << endl; 
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
}

extern "C" void IDL_get_intensity_autocorrelation(int argc, void * argv[]){
  check_objects();
  int nx = esw->get_size_x();
  int ny = esw->get_size_y();
  
  Double_2D temp(nx,ny);
  reco->get_intensity_autocorrelation(temp);
  copy_from_double_2d(temp,(double*) argv[0]);
}

extern "C" void IDL_get_support(int argc, void * argv[]){
  check_objects();
  int nx = esw->get_size_x();
  int ny = esw->get_size_y();
  
  Double_2D temp(nx,ny);
  reco->get_support(temp);
  copy_from_double_2d(temp,(double*) argv[0]); 

}


extern "C" void IDL_get_error(int argc, void * argv[]){
  check_objects();
  *(double*) argv[0] = reco->get_error();
}

extern "C" void IDL_get_transmission_function(int argc, void * argv[]){
  check_objects();
  if(typeid(*reco)!=typeid(FresnelCDI)){

    ostringstream oss (ostringstream::out);
    oss << "Sorry, can't get the transmission function for "
	<< "anything other than "<< typeid(FresnelCDI).name() <<" reconstuction. "
	<< "You are doing "<<typeid(*reco).name() 
	<< " reconstruction." << endl;
    IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());    
    return;
  }
  
  int nx = esw->get_size_x();
  int ny = esw->get_size_y();

  Complex_2D temp(nx,ny);
  ((FresnelCDI*) reco)->get_transmission_function(temp);
  copy_from_complex_2d(temp,(IDL_COMPLEX*) argv[0]); 
}


extern "C" void IDL_print_algorithm(int argc, void * argv[]){
  check_objects();

  //now we do something tricky to redirect stdout into a string
  //this is necessary because IDL doesn't play well with stdout in
  //windows.
  streambuf * backup = cout.rdbuf(); //store a pointer to the stdout buffer
  streambuf * str_buffer = new stringbuf(); //create a new string buffer
  cout.rdbuf(str_buffer);  //redirect stdout to the string buffer

  reco->print_algorithm(); //get the algorithm

  cout.rdbuf(backup);  //restore cout to stdout

  //pass the string to IDL
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, ( (stringbuf*) str_buffer)->str().c_str());

  //clean up
  delete str_buffer;
}


extern "C" void IDL_propagate_from_detector(int argc, void * argv[]){
  check_objects();
  Complex_2D temp(esw->get_size_x(),esw->get_size_y());
  copy_to_complex_2d(temp,(IDL_COMPLEX*) argv[0]); 
  reco->propagate_from_detector(temp);
  copy_from_complex_2d(temp,(IDL_COMPLEX*) argv[1]); 
}


extern "C" void IDL_propagate_to_detector(int argc, void * argv[]){
  check_objects();
  Complex_2D temp(esw->get_size_x(),esw->get_size_y());
  copy_to_complex_2d(temp,(IDL_COMPLEX*) argv[0]); 
  reco->propagate_to_detector(temp);
  copy_from_complex_2d(temp,(IDL_COMPLEX*) argv[1]); 
}

extern "C" void IDL_apply_support(int argc, void * argv[]){
  check_objects();
  Complex_2D temp(esw->get_size_x(),esw->get_size_y());
  copy_to_complex_2d(temp,(IDL_COMPLEX*) argv[0]); 
  reco->apply_support(temp);
  copy_from_complex_2d(temp,(IDL_COMPLEX*) argv[1]); 
}

extern "C" void IDL_scale_intensity(int argc, void * argv[]){
  check_objects();
  Complex_2D temp(esw->get_size_x(),esw->get_size_y());
  copy_to_complex_2d(temp,(IDL_COMPLEX*) argv[0]); 
  reco->scale_intensity(temp);
  copy_from_complex_2d(temp,(IDL_COMPLEX*) argv[1]); 
}
