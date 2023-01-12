// Copyright 2011 Nadia Davidson 
// for The ARC Centre of Excellence in Coherent X-ray Science. 
//
// This program is distributed under the GNU General Public License. 
// We also ask that you cite this software in publications where you made 
// use of it for any part of the data analysis.

//standard headers
#include <cstdlib>
#include <iostream>
#include <typeinfo>
#include <sstream>
#include <string>

//idl header
#include "idl_export.h"

//NADIA software headers
#include <Complex_2D.h>
#include <Double_2D.h>
#include <BaseCDI.h>
#include <PlanarCDI.h>
#include <FresnelCDI.h>
#include <FresnelCDI_WF.h>
#include <PartialCDI.h>
#include <PartialCharCDI.h>
#include <PolyCDI.h>
#include <PhaseDiverseCDI.h>
#include <TransmissionConstraint.h>
#include <io.h>

Complex_2D* esw = 0;
BaseCDI * reco = 0;
TransmissionConstraint * constraint = 0;

PhaseDiverseCDI * phaseDiv = 0;
std::vector<Complex_2D*> phaseDiv_complex;
std::vector<BaseCDI*> phaseDiv_cdi;
std::vector<TransmissionConstraint*> phaseDiv_constraint;

int total_iters = 0;
using namespace std;

/************ helper methods ******************/

//methods for converting between the IDL types and the C++ types.
void copy_to_double_2d(Double_2D & nadia_array, double * IDL_array){

  int ny = nadia_array.get_size_x();
  int nx = nadia_array.get_size_y();

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      nadia_array.set(j,nx-1-i,*(IDL_array));
      IDL_array++;
    }
  }
}

void copy_from_double_2d(const Double_2D & nadia_array, double * IDL_array){

  int ny = nadia_array.get_size_x();
  int nx = nadia_array.get_size_y();

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      *(IDL_array)= nadia_array.get(j,nx-1-i);
      IDL_array++;
    }
  }
}

void copy_to_complex_2d(Complex_2D & nadia_array, IDL_COMPLEX * IDL_array){

  int ny = nadia_array.get_size_x();
  int nx = nadia_array.get_size_y();

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      nadia_array.set_real(j,nx-1-i,(*IDL_array).r);
      nadia_array.set_imag(j,nx-1-i,(*IDL_array).i);
      IDL_array++;
    }
  }
}

void copy_from_complex_2d(const Complex_2D & nadia_array, IDL_COMPLEX * IDL_array){

  int ny = nadia_array.get_size_x();
  int nx = nadia_array.get_size_y();  

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      (*IDL_array).r= nadia_array.get_real(j,nx-1-i);
      (*IDL_array).i= nadia_array.get_imag(j,nx-1-i);
      IDL_array++;
    }
  }
}

//error checking method. Make sure the objects have memory allocated.
void check_objects(){
  if(!esw || !reco){
    char buf[] = "You need to call one of the NADIA_INIT_ routines before calling this routine";
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
  int status = read_dbin(filename.s,nx,ny,temp);
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


/************ reconstruction methods ******************/



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

  if(constraint!=0){
    constraint->delete_complex_constraint_regions();
    delete constraint;
    constraint = 0;
  }

  if(phaseDiv!=0){
    delete phaseDiv;
    phaseDiv = 0;
  }

  while(!phaseDiv_complex.empty()){
    delete phaseDiv_complex.back();
    phaseDiv_complex.pop_back();
  }

  while(!phaseDiv_cdi.empty()){
    delete phaseDiv_cdi.back();
    phaseDiv_cdi.pop_back();
  }

  while(!phaseDiv_constraint.empty()){
    phaseDiv_constraint.back()->delete_complex_constraint_regions();
    delete phaseDiv_constraint.back();
    phaseDiv_constraint.pop_back();
  }

}


void common_init(int argc, void * argv[], int max_args){

  if(phaseDiv==0)
    IDL_deallocate_memory(0,0);

  total_iters = 0;

  int ny = *(int*) argv[0];
  int nx = *(int*) argv[1];

//  std::cout<<"nx is "<<nx<<" and ny is "<<ny<<"\n\n";

  esw = new Complex_2D(nx,ny);  

  if(argc == max_args){
    copy_to_complex_2d(*esw,(IDL_COMPLEX*) argv[max_args-1]);
  }

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

extern "C" void IDL_partial_init(int argc, void * argv[])
{

  common_init(argc, argv, 11);

  reco = new PartialCDI(*esw,
      *(double*) argv[2],
      *(double*) argv[3],
      *(double*) argv[4],
      *(double*) argv[5],
      *(double*) argv[6],
      *(double*) argv[7],
      *(double*) argv[8],
      *(double*) argv[9]);

}

extern "C" void IDL_part_char_init(int argc, void * argv[])
{

  common_init(argc, argv, 8);

  reco = new PartialCharCDI(*esw,
      *(double*) argv[2],
      *(double*) argv[3],
      *(double*) argv[4],
      *(double*) argv[5],
      *(double*) argv[6]);

}

extern "C" void IDL_poly_init(int argc, void * argv[])
{

  common_init(argc, argv, 3);
  reco = new PolyCDI(*esw);

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

//this doesn't actually call any of the NADIA library code
extern "C" void IDL_get_round_support(int argc, void * argv[]){

  int ny = *(int*) argv[0];
  int nx = *(int*) argv[1];

  double radius = *(double*) argv[2]; 
  Double_2D temp_2D(nx,ny);
  double i0 = (nx-1)/2.0;
  double j0 = (ny-1)/2.0;
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      if(sqrt((i-i0)*(i-i0)+(j-j0)*(j-j0)) 
	  < radius*sqrt(j0*j0+i0*i0))
	temp_2D.set(i,j,100);
    }
  }
  copy_from_double_2d(temp_2D,(double*) argv[3]); 

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
  int alg = BaseCDI::getAlgFromName(alg_name.s);
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

  if(typeid(*reco)!=typeid(PlanarCDI)){

    ostringstream oss (ostringstream::out);
    oss << "Sorry, can't get the autocorrelation function for "
      << "anything other than "<< typeid(PlanarCDI).name() <<" reconstruction. "
      << "You are doing "<<typeid(*reco).name() 
      << " reconstruction." << endl;
    IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());    
    return;
  }

  int nx = esw->get_size_x();
  int ny = esw->get_size_y();

  Double_2D temp(nx,ny);
  ((PlanarCDI*)reco)->get_intensity_autocorrelation(temp);
  copy_from_double_2d(temp,(double*) argv[0]);
}

extern "C" void IDL_get_support(int argc, void * argv[]){
  check_objects();
  int nx = esw->get_size_x();
  int ny = esw->get_size_y();

  //  Double_2D temp(nx,ny);
  //reco->get_support(temp);
  copy_from_double_2d(reco->get_support(),(double*) argv[0]); 

}


extern "C" void IDL_get_error(int argc, void * argv[]){
  check_objects();
  *(double*) argv[0] = reco->get_error();
}

extern "C" void IDL_get_transmission_function(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(FresnelCDI)){

    int nx = esw->get_size_x();
    int ny = esw->get_size_y();

    Complex_2D temp(nx,ny);
    ((FresnelCDI*) reco)->get_transmission_function(temp);
    copy_from_complex_2d(temp,(IDL_COMPLEX*) argv[0]); 
    return;
  }

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't get the transmission function for "
    << "anything other than "<< typeid(PartialCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return;
}

extern "C" void IDL_get_transmission(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(PartialCDI)){

    int nx = esw->get_size_x();
    int ny = esw->get_size_y();

    Complex_2D temp(nx,ny);
    temp=((PartialCDI*) reco)->get_transmission();
    copy_from_complex_2d(temp,(IDL_COMPLEX*) argv[0]);
    return;
  }

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't get the transmission function for "
    << "anything other than "<< typeid(PartialCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return;
}

extern "C" void IDL_set_transmission(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(PartialCDI)){

    int ny = *(int*) argv[0];
    int nx = *(int*) argv[1];

    Complex_2D temp_2D(nx,ny);

    copy_to_complex_2d(temp_2D, (IDL_COMPLEX*) argv[2]);

    ((PartialCDI*) reco)->set_transmission(temp_2D);
    return;
  }

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't set the transmission function for "
    << "anything other than "<< typeid(PartialCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return;
}

extern "C" void IDL_initialise_matrices(int argc, void * argv[]){

  //check_objects();
  double dnleg = *(double*) argv[0];
  double dnmodes = *(double*) argv[1];

  if(dnleg < dnmodes){
    ostringstream oss (ostringstream::out);
    oss <<"The order of the legendre approximation"
      <<" must be greater than or equal to the number"
      <<" of modes. The matrices have not been"
      <<" initialised" <<endl;

    IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
    return;
  }

  if(typeid(*reco)==typeid(PartialCDI)){

    int nleg = (int) dnleg;
    int nmodes = (int) dnmodes;

    ((PartialCDI*) reco)->initialise_matrices(nleg, nmodes);

  }else{

    ostringstream oss (ostringstream::out);
    oss << "Sorry, can't initialise the matrices for "
      << "anything other than "<< typeid(PartialCDI).name() <<" reconstruction. "
      << "You are doing "<<typeid(*reco).name()
      << " reconstruction." << endl;
    IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
    return;
  }
}

extern "C" void IDL_set_spectrum(int argc, void * argv[]){

  if(typeid(*reco)==typeid(PolyCDI)){
    IDL_STRING filename = *(IDL_STRING*)argv[0];
    ((PolyCDI*) reco)->set_spectrum(filename.s);
  }else{
    ostringstream oss (ostringstream::out);
    oss << "Sorry, can't set the spectrum for "
      << "anything other than "<< typeid(PolyCDI).name() <<" reconstruction. "
      << "You are doing "<<typeid(*reco).name()
      << " reconstruction." << endl;
    IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
    return;
  }
}

extern "C" void IDL_get_mode(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(PartialCDI)){

    int nx = esw->get_size_x();
    int ny = esw->get_size_y();

    double dnmode = *(double*)argv[1];

    int nmode = (int) dnmode;


    Complex_2D temp(nx,ny);

    temp=((PartialCDI*) reco)->get_mode(nmode);
    copy_from_complex_2d(temp,(IDL_COMPLEX*) argv[0]);
    return;
  }

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't get the modes function for "
    << "anything other than "<< typeid(PartialCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return;
}

extern "C" void IDL_set_threshold(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(PartialCDI)){

    double thresh = *(double*) argv[0];

    ((PartialCDI*) reco)->set_threshold(thresh);
    return;
  }

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't set the threshhold for "
    << "anything other than "<< typeid(PartialCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return;
}

extern "C" void IDL_set_initial_coherence_guess(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(PartialCharCDI)){

    double lx = *(double*) argv[0];
    double ly = *(double*) argv[1];

    ((PartialCharCDI*) reco)->set_initial_coherence_guess(lx, ly);
    return;
  }

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't set the initial coherence guess for "
    << "anything other than "<< typeid(PartialCharCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return;
}

extern "C" void IDL_set_initial_coherence_guess_in_m(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(PartialCharCDI)){

    double lx = *(double*) argv[0];
    double ly = *(double*) argv[1];

    ((PartialCharCDI*) reco)->set_initial_coherence_guess_in_m(lx, ly);
    return;           
  }

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't set the initial coherence guess for "
    << "anything other than "<< typeid(PartialCharCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return;     
}

extern "C" void IDL_set_minima_search_bounds_coefficient(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(PartialCharCDI)){

    double coef = *(double*) argv[0];

    ((PartialCharCDI*) reco)->set_minima_search_bounds_coefficient(coef);
    return;
  }

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't set the minima search bound coefficient for "
    << "anything other than "<< typeid(PartialCharCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return;
}

extern "C" void IDL_set_minima_search_tolerance(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(PartialCharCDI)){

    double tol = *(double*) argv[0];

    ((PartialCharCDI*) reco)->set_minima_search_tolerance(tol);
    return;             
  }                       

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't set the minima search tolerance for "
    << "anything other than "<< typeid(PartialCharCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return;     
} 

extern "C" void IDL_set_minima_search_tolerance_in_m(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(PartialCharCDI)){

    double tol = *(double*) argv[0];

    ((PartialCharCDI*) reco)->set_minima_search_tolerance_in_m(tol);
    return;            
  }

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't set the minima search tolerance for "
    << "anything other than "<< typeid(PartialCharCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return;     
}

extern "C" void IDL_set_minima_moving_average_weight(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(PartialCharCDI)){

    double w = *(double*) argv[0];

    ((PartialCharCDI*) reco)->set_minima_moving_average_weight(w);
    return;
  }

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't set the minima average weight for "
    << "anything other than "<< typeid(PartialCharCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return;
}

extern "C" void IDL_set_minima_recalculation_interval(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(PartialCharCDI)){

    int ival = *(int*) argv[0];

    ((PartialCharCDI*) reco)->set_minima_recalculation_interval(ival);
    return;
  }

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't set the minima recalculation interval "
    << "anything other than "<< typeid(PartialCharCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return;
}

extern "C" IDL_LONG IDL_get_x_coherence_length(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(PartialCharCDI)){

    return ((PartialCharCDI*) reco)->get_x_coherence_length();

  }

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't get the coherence length for "
    << "anything other than "<< typeid(PartialCharCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return 0;
}

extern "C" IDL_LONG IDL_get_y_coherence_length(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(PartialCharCDI)){

    return ((PartialCharCDI*) reco)->get_y_coherence_length();
  }

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't get the coherence length for "
    << "anything other than "<< typeid(PartialCharCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return 0;
}

extern "C" IDL_LONG IDL_get_x_coherence_length_in_pixels(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(PartialCharCDI)){

    return ((PartialCharCDI*) reco)->get_x_coherence_length_in_pixels();

  }

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't get the coherence length for "
    << "anything other than "<< typeid(PartialCharCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return 0;
}

extern "C" IDL_LONG IDL_get_y_coherence_length_in_pixels(int argc, void * argv[]){
  check_objects();

  if(typeid(*reco)==typeid(PartialCharCDI)){

    return ((PartialCharCDI*) reco)->get_y_coherence_length_in_pixels();

  }               

  ostringstream oss (ostringstream::out);
  oss << "Sorry, can't get the coherence length for "
    << "anything other than "<< typeid(PartialCharCDI).name() <<" reconstruction. "
    << "You are doing "<<typeid(*reco).name()
    << " reconstruction." << endl;
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, oss.str().c_str());
  return 0;     
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



//------------------------------------------------------------//
// Complex contraint code
//------------------------------------------------------------//

//args: true/false
extern "C" void IDL_set_trans_unity_constraint(int argc, void * argv[]){
  check_objects();
  if(!constraint){
    constraint = new TransmissionConstraint();
    reco->set_complex_constraint(*constraint);
  }

  constraint->set_enforce_unity( *(bool*) argv[0] );

}

//args: true/false 
extern "C" void IDL_set_charge_flipping(int argc, void * argv[]){
  check_objects();
  if(!constraint){
    constraint = new TransmissionConstraint();
    reco->set_complex_constraint(*constraint);
  }

  constraint->set_charge_flipping( *(bool*) argv[0] );

}

//args: region, alpha1, alpha2, fixed c (optional)
extern "C" void IDL_add_complex_constraint_region(int argc, void * argv[]){

  //first check that the appropriate memory has already been allocated
  check_objects();
  if(!constraint){
    constraint = new TransmissionConstraint();
    reco->set_complex_constraint(*constraint);
  }

  //copy the input into the right format.
  Double_2D temp_region(reco->get_size_x(),reco->get_size_y());
  copy_to_double_2d(temp_region,(double*) argv[0]); 
  double alpha1 = *(double*) argv[1];
  double alpha2 = *(double*) argv[2];

  //create the new constraint
  ComplexConstraint * new_region;
  new_region = new ComplexConstraint(temp_region, alpha1, alpha2);
  constraint->add_complex_constraint(*new_region);
  //  constraint_regions.push_back(new_region);

  //if a fixed value for c=beta/delta has been given:
  if(argc==4)
    new_region->set_fixed_c(*(double*) argv[3]);

}



//------------------------------------------------------------//
// Binding to the PhaseDiverseCDI functions
//------------------------------------------------------------//

void check_pd(){
  if(phaseDiv==0){
    IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO,
	"ERROR: trying to use phase diverse functions before initialising.");
    exit(0);
  }
}

void check_pd_trans(){

  check_pd();

  if(phaseDiv && phaseDiv->get_transmission()==0){
    IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO,
	"ERROR: trying to use phase diverse functions before adding frames.");
    exit(0);
  }
}


extern "C" void IDL_phase_diverse_init(int argc, void * argv[]){

  double beta = *(double*) argv[0];
  double gamma = *(double*) argv[1];
  bool parallel = *(bool*) argv[2];

  if(phaseDiv==0)
    phaseDiv = new PhaseDiverseCDI(beta, gamma, parallel);
  else
    IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO,
	"Warning: trying to intialise the phase diverse code more than once.");

}

//double x=0, double y=0, double alpha=1
extern "C" void IDL_phase_diverse_add_position(int argc, void * argv[]){
  check_pd();

  double x = *(double*) argv[0];
  double y = *(double*) argv[1];
  double alpha = *(double*) argv[2];

  BaseCDI * local = reco;
  Complex_2D * local_esw = esw;
  TransmissionConstraint * local_con = constraint;

  //redirect stardard out to a buffer and put in IDL message format
  streambuf * backup = cout.rdbuf(); //store a pointer to the stdout buffer
  streambuf * str_buffer = new stringbuf(); //create a new string buffer
  cout.rdbuf(str_buffer);  //redirect stdout to the string buffer

  //add the new position
  phaseDiv->add_new_position(local, x, y, alpha);

  cout.rdbuf(backup);  //restore cout to stdout
  //pass the string to IDL
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, ( (stringbuf*) str_buffer)->str().c_str());
  //clean up
  delete str_buffer;


  //book keeping..
  phaseDiv_cdi.push_back(local);
  phaseDiv_complex.push_back(local_esw);
  phaseDiv_constraint.push_back(local_con);

  reco = 0;
  esw = 0;
  constraint = 0;

}

extern "C" void IDL_phase_diverse_init_estimate(int argc, void * argv[]){

  check_pd_trans();
  phaseDiv->initialise_estimate();

  //return the result
  Complex_2D * result = phaseDiv->get_transmission();
  copy_from_complex_2d(*result ,(IDL_COMPLEX*) argv[0]); 

}

extern "C" void IDL_phase_diverse_iterate(int argc, void * argv[]){
  check_pd_trans();
  int iterations = *(IDL_LONG*) argv[0];

  for(int i=0; i<iterations; i++){

    //redirect stardard out to a buffer and put in IDL message format
    streambuf * backup = cout.rdbuf(); //store a pointer to the stdout buffer
    streambuf * str_buffer = new stringbuf(); //create a new string buffer
    cout.rdbuf(str_buffer);  //redirect stdout to the string buffer

    //do the actual iteration.
    phaseDiv->iterate(); 

    cout.rdbuf(backup);  //restore cout to stdout
    //pass the string to IDL
    IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, ( (stringbuf*) str_buffer)->str().c_str());
    //clean up
    delete str_buffer;
  }

  //return the result
  Complex_2D * result = phaseDiv->get_transmission();
  copy_from_complex_2d(*result ,(IDL_COMPLEX*) argv[1]); 
}

//Allows the dimensions of the image to be passed back to the IDL code
extern "C" IDL_LONG IDL_get_phase_diverse_array_x_size(int argc, void *argv[]){
  check_pd_trans();
  return phaseDiv->get_transmission()->get_size_x();
}

extern "C" IDL_LONG IDL_get_phase_diverse_array_y_size(int argc, void *argv[]){
  check_pd_trans();
  return phaseDiv->get_transmission()->get_size_y();
}


extern "C" void IDL_phase_diverse_iterations_per_cycle(int argc, void * argv[]){
  check_pd();
  int iterations = *(IDL_LONG*) argv[0];
  phaseDiv->set_iterations_per_cycle(iterations);
}


extern "C" void IDL_phase_diverse_set_transmission(int argc, void * argv[]){
  check_pd();
  int nx = *(IDL_LONG*) argv[0];
  int ny = *(IDL_LONG*) argv[1];

  Complex_2D new_transmission(nx,ny);
  copy_to_complex_2d(new_transmission ,(IDL_COMPLEX*) argv[2]); 

  phaseDiv->set_transmission(new_transmission);
}

extern "C" void IDL_phase_diverse_adjust_positions(int argc, void * argv[]){
  check_pd_trans();

  int type = *(IDL_LONG*) argv[0]; //0 - cross correlation, 1 - error minimisation.
  bool forward = *(bool*) argv[1];
  int x_min = *(IDL_LONG*) argv[2]; 
  int x_max = *(IDL_LONG*) argv[3]; 
  int y_min = *(IDL_LONG*) argv[4]; 
  int y_max = *(IDL_LONG*) argv[5]; 
  double step_size = *(double*) argv[6];

  //redirect stardard out to a buffer and put in IDL message format
  streambuf * backup = cout.rdbuf(); //store a pointer to the stdout buffer
  streambuf * str_buffer = new stringbuf(); //create a new string buffer
  cout.rdbuf(str_buffer);  //redirect stdout to the string buffer


  //adjust the positions
  phaseDiv->adjust_positions(type,forward, 
      x_min, x_max, 
      y_min, y_max, 
      step_size);


  cout.rdbuf(backup);  //restore cout to stdout
  //pass the string to IDL
  IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, ( (stringbuf*) str_buffer)->str().c_str());
  //clean up
  delete str_buffer;


}

extern "C" double IDL_phase_diverse_get_final_x_position(int argc, void * argv[]){
  check_pd_trans();
  //*(double*) argv[1] = phaseDiv->get_final_x_position (*(IDL_LONG*)argv[0] );
  return phaseDiv->get_final_x_position (*(IDL_LONG*)argv[0] );
}

extern "C" double IDL_phase_diverse_get_final_y_position(int argc, void * argv[]){
  check_pd_trans();
  return phaseDiv->get_final_y_position( *(IDL_LONG*)argv[0] );
}



