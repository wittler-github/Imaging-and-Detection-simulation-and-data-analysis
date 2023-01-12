// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <cstdlib> 
#include <Complex_2D.h>
#include <Double_2D.h>
#include <FresnelCDI.h>
#include <PlanarCDI.h>
#include <PhaseDiverseCDI.h>
#include <io.h>
#include <iomanip>
#include <sstream>
#include <typeinfo>
#include <utils.h>

using namespace std;


//constructor for the class which handles phase diversity.
PhaseDiverseCDI::PhaseDiverseCDI(
				 double beta, 
				 double gamma, 
				 bool parallel,
				 int granularity):scale(granularity),
						  beta(beta), 
						  gamma(gamma), 
						  parallel(parallel),
						  nx(0),
						  ny(0),
						  x_min(0),
						  y_min(0),
						  weights_set(false),
                                                  iterations_per_cycle(1),
                                                  object(0){
  
};


//destructor for cleaning up
PhaseDiverseCDI::~PhaseDiverseCDI(){

  while(single_result.size()!=0){
    Complex_2D * temp = single_result.back();
    single_result.pop_back();
    delete temp;
    
    Double_2D * temp_d = weights.back();
    weights.pop_back();
    delete temp_d;
  }
  
  if(object)
    delete object;
  
}

//change the size of the global transmission function 
//this function is used for dynamically determining the size
//of the glocal transmission function array.
void PhaseDiverseCDI::reallocate_object_memory(int new_nx,int new_ny){

  if(object)
    delete object;

  object = new Complex_2D(new_nx,new_ny);
  
  nx = new_nx;
  ny = new_ny;

}

//return the global transmision function.
Complex_2D * PhaseDiverseCDI::get_transmission(){
  return object;
}


//set the global transmission function (note that
//local transmission function will not be automatically updated).
//this function can be used to set the initial estimate.
void PhaseDiverseCDI::set_transmission(Complex_2D & new_transmission){

  int new_nx = new_transmission.get_size_x();
  int new_ny = new_transmission.get_size_y();

  if(!object)
    reallocate_object_memory(new_nx,new_ny);

  if(new_nx!=nx||new_ny!=ny){
    cout << "Can not set the transmission function in "
	 << "PhaseDiverseCDI because the complex array "
	 << "given does not have the same dimensions as "
	 << "the current transmission function"<<endl;
    return;
  }

  object->copy(new_transmission);
  
}

//add a new frame to the reconstruction
//the size of the global sample will be adjusted automatically.
void PhaseDiverseCDI::add_new_position(BaseCDI * local,
				       double x, 
				       double y,
				       double alpha){
  
  weights_set = false;

  //get the size of the frame which is being added
  int lnx = local->get_size_x() ;
  int lny = local->get_size_y() ;
  
  //add it to the list of frames with it's position.
  singleCDI.push_back(local);
  x_position.push_back(x);
  y_position.push_back(y);
  this->alpha.push_back(alpha);
  single_result.push_back(new Complex_2D(lnx,lny));
  weights.push_back(new Double_2D(lnx,lny));

  cout << "Added position "<<singleCDI.size()-1<<endl;
  
  //if this is the first frame we will need to create the 
  //global transmission object as well
  if(!object){
    x_min = -x;
    y_min = -y;
    reallocate_object_memory(lnx*scale,lny*scale);
  }

  //dynamically increase the object size to fit this frame in

  //what is the position (globally) of the first 
  //pixel if we don't increase the frame size  
  int extra_x=0;
  int extra_y=0;

  int global_x_min = get_global_x_pos(0,x)*scale;  
  int global_y_min = get_global_y_pos(0,y)*scale;  
  int global_x_max = get_global_x_pos(lnx,x)*scale; 
  int global_y_max = get_global_y_pos(lny,y)*scale; 

  if(global_x_min<0){
    x_min = -x;
    extra_x += -global_x_min;
  }
  if(global_y_min<0){
    y_min = -y;
    extra_y += -global_y_min;
  }
  if(global_x_max>nx)
    extra_x += global_x_max-nx;
  if(global_y_max>ny)
    extra_y += global_y_max-ny;

  //if required, increase the global object size
  if(extra_x || extra_y)
    reallocate_object_memory(nx+extra_x, ny+extra_y);

}



//initialise the global estimate based on the starting point
//for each sub-frame.
void PhaseDiverseCDI::initialise_estimate(){

  //start by setting the magnitude to 1 and phase to 0
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      object->set_real(i,j,1);
      object->set_imag(i,j,0);
    }
  }
  

  //add the transmission from each frame
  //to the global object.
  for(int i=0; i<singleCDI.size(); i++){
    add_to_object(i);
  }

 
}

//do a 'big' iteration.
/*void*/double * PhaseDiverseCDI::iterate(){

  double error_sum_c = 0; //For article
  double error_sum_m = 0;//article
  double error_sum_p = 0;//article
  double mean_error_c = 0;//For article
  double mean_error_m = 0;//For article
  double mean_error_p = 0;//For article
  //double max_iterations = 10; //For article
  
  static int total_iterations = 0;

  cout << "Iteration "<<total_iterations<<endl;
  cout << "Iterations per cycle " << iterations_per_cycle << endl; //for article

  //for each frame do a 'small' iteration.
  for(int i=0; i<singleCDI.size(); i++){

    int x, y;

    //if(i!=0 && (total_iterations==2 || total_iterations==4))
    //  check_position(i);
    update_from_object(i);

    //if the user wanted do a few iterations before updating to the
    //global object.
    for(int j=0; j<iterations_per_cycle; j++){
				     
      //temp
      /**      static bool flag = false;
      if(!flag){
	char buf[50];
	Double_2D result(nx,ny);
	Complex_2D trans(nx,ny);
	sprintf(buf,"itr_%i_%i_1.tiff",i,j);
	//      ((FresnelCDI*)singleCDI.at(i))->get_transmission_function(trans);
	trans.copy(*get_transmission());
	trans.get_2d(MAG,result);
	write_image(buf,result,false,0,1.0);
	flag=true;
	}**/
      
      singleCDI.at(i)->iterate();
      
      cout << "Error for frame "<<i<<" for iteration " << j << " is "
	   << singleCDI.at(i)->get_error() << endl;

     if( j == iterations_per_cycle-1 ){ 
     error_sum_m += singleCDI.at(i)->get_error();
     error_sum_p += singleCDI.at(i)->get_error_p();
     error_sum_c += singleCDI.at(i)->get_error_c();
				     }

    }

    //if we are running in series mode, update the
    //small transmission to the large transmission function
    //after each 'small' iteration
    if(!parallel)
      add_to_object(i);

  }
  
  //if we are running in parallel mode, wait until at least one
  //'small' iteration have been performed for each frame. Then
  //merge the results to form the new global transmission function. 
  if(parallel){

    scale_object(1-beta);

    for(int i=0; i<singleCDI.size(); i++){
      add_to_object(i);
    }
  }

  total_iterations++;

  //if(max_iterations == total_iterations){total_iterations = 0;} //For article
  mean_error_c = error_sum_c/(singleCDI.size()); //For article
  mean_error_m = error_sum_m/(singleCDI.size()); //For article
  mean_error_p = error_sum_p/(singleCDI.size()); //For article
  cout << "Mean Chi Squared " << mean_error_c << " " << mean_error_m << " " << mean_error_p << endl;
  
  double * mean_error_cmp_array = new double[3];//article
  /*double  mean_error_cmp[3];
         mean_error_cmp_array = mean_error_cmp;*/
	 mean_error_cmp_array[0] = mean_error_c;
	 mean_error_cmp_array[1] = mean_error_m;
	 mean_error_cmp_array[2] = mean_error_p;
	

  return mean_error_cmp_array; //For article
 
   
};


//scale each element in the global transmission object.
//note that the region outside the sample will not be scaled.
//this function is only used when performing the reconstruction
//in parallel mode.
void PhaseDiverseCDI::scale_object(double factor){

  set_up_weights();

  int frames = singleCDI.size();

  //loop over each element in the global transmission function.
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){

      bool not_in_image = true;
      int n_probe = 0;

      //loop over each frame and check their weight function to determine
      //if the current pixel is included inside the sample or not.
      while(not_in_image&&n_probe<frames){

	int i_ = get_local_x_pos(i/scale, x_position.at(n_probe));  
	int j_ = get_local_y_pos(j/scale, y_position.at(n_probe));
	
	if(i_>=0&&j_>=0
	   &&i_<weights.at(n_probe)->get_size_x()
	   &&j_<weights.at(n_probe)->get_size_y())
	  not_in_image=(weights.at(n_probe)->get(i_,j_)==0); 
	
	n_probe++;	
      }

      //in the case the pixel is outside, set the magntiude 
      //to 1 and phase to 0.
      if(not_in_image){
	object->set_real(i,j,1.0);
	object->set_imag(i,j,0.0);
      } //otherwise scale the value of the pixel.
      else{
	object->set_real(i,j,factor*object->get_real(i,j));
	object->set_imag(i,j,factor*object->get_imag(i,j));
      }	
      
    }
  }
  
}

//update the array 'result' with the contents of either the
//transmission function (for FresnelCDI) of the exit-surface-wave
//(for plane-wave).
void PhaseDiverseCDI::get_result(BaseCDI * local, Complex_2D & result){
  if(typeid(*local)==typeid(FresnelCDI)){
    ((FresnelCDI*)local)->get_transmission_function(result);
  }
  else
    result.copy(local->get_exit_surface_wave());
};

//copy the array 'result' back to the CDI object 
void PhaseDiverseCDI::set_result(BaseCDI * local, Complex_2D & result){
  
  if(typeid(*local)==typeid(FresnelCDI)){
    ((FresnelCDI*)local)->set_transmission_function(result);
  }
  else
    local->set_exit_surface_wave(result);
  
}

//frames need to be in order for this to work.
void PhaseDiverseCDI::adjust_positions(int type, bool forward, 
				       int min_x, int max_x,
				       int min_y, int max_y,
				       double step_size){

  if(type!=CROSS_CORRELATION && type!=MINIMUM_ERROR){
    cerr << "The alignment type given to"
	 << " PhaseDiverseCDI::adjust_positions is not"
	 << " a known type. Please consider using"
	 << " PhaseDiverseCDI::CROSS_CORRELATION or "
	 << " PhaseDiverseCDI::MINIMUM_ERROR." 
	 << endl;
    return;
  }

  bool parallel_c = parallel; 

  //if we are running in series, switch to 
  //parallel mode and do 1 full iteration:
  /** if(!parallel){
    parallel = true;
    weights_set=false;
    for(int i=0; i<10; i++)
      iterate();
      }**/

  //make a copy of the tranmission function
  Complex_2D * pointer_to_object = object;
  int nx_c = nx;
  int ny_c = ny;
  double x_min_c = x_min; 
  double y_min_c = y_min; 

  double beta_c = beta;

  //reset the parameters to make it work in parrallel with weights of one.
  parallel = false;

  //leave the first frame, as everything will be aligned to it.
  int n=1;
  int limit=singleCDI.size();
  if(!forward){
    n=singleCDI.size()-2;
    limit=-1;
  }

  Complex_2D temp_object(nx,ny);
  object = &temp_object;
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      temp_object.set_real(i,j,1);
      temp_object.set_imag(i,j,0);
    }
  }

  beta = 1.0;
  weights_set=false;

  double n_first_frame = 0;

  if(!forward)
    n_first_frame = singleCDI.size()-1;

  singleCDI.at(n_first_frame)->iterate();
  get_result(singleCDI.at(n_first_frame),*(single_result.at(n_first_frame)));
  add_to_object(n_first_frame);
  
  while(n!=limit){

    //store some information for later use
    double before_x = x_position.at(n);
    double before_y = y_position.at(n);

    int lnx = single_result.at(0)->get_size_x();
    int lny = single_result.at(0)->get_size_y();

    //do an iteration of the frame on it's own.
    singleCDI.at(n)->iterate();
    get_result(singleCDI.at(n),*(single_result.at(n)));

    if(type==CROSS_CORRELATION){

      //this array will store a copy of the local frame with 1-mag.
      Double_2D temp_single(lnx,lny);

      //this array will store a copy of the global frame with 1-mag.
      //in a sub-grid.
      Double_2D temp_others(lnx*scale,lny*scale);
      
      for(int i=0; i<lnx*scale; i++){
	for(int j=0; j<lny*scale; j++){

	  int i_ = i - before_x*scale -  x_min*scale;
	  int j_ = j - before_y*scale -  y_min*scale;

	  if( i_>=0 && j_>=0 && i_<nx && j_<ny ){
	    double obj_mag = object->get_mag(i_,j_);
	    if(obj_mag<1)
	      temp_others.set(i,j,1-obj_mag);
	  }
	  if(i<lnx && j<lny && single_result.at(n)->get_mag(i,j) <1 ){
	    double new_value = 1-single_result.at(n)->get_mag(i,j);
	    temp_single.set(i,j,new_value);
	  }
	}
      }
      

      /**      char name[80];
      Double_2D pic(nx,ny);
      object->get_2d(MAG,pic);
      sprintf(name,"pic_%i.tiff",n);
      write_image(name,pic,false, 0,1);
      **/

      Double_2D * temp_single_int = 0; 
      Double_2D * weight_single_int = new Double_2D(lnx*scale,lny*scale);
      
      if(scale==1){
	temp_single_int = &temp_single;
	weight_single_int->copy(singleCDI.at(n)->get_support());
      }
      else{
	temp_single_int = new Double_2D(lnx*scale,lny*scale);
	interpolate(temp_single,*temp_single_int);
	interpolate(singleCDI.at(n)->get_support(),*weight_single_int);
      }

      /**      static int counter = 0;
      char buff[90];
      sprintf(buff,"temp_single_%i.tiff",counter);
      write_image(buff,*temp_single_int,false);
      //cout << "not here"<<endl;
      
      sprintf(buff,"temp_object_%i.tiff",counter);
      write_image(buff,temp_others,false);
      counter++; **/

      int new_x=0;
      int new_y=0;

      align(*temp_single_int, temp_others, 
	    new_x, new_y,
	    min_x, max_x,
	    min_y, max_y,
	    weight_single_int,
	    &temp_others);
      
      x_position.at(n) = before_x + new_x/((double)scale);
      y_position.at(n) = before_y + new_y/((double)scale);

      /**      cout << "moving prob " << n << " from ("
	   << before_x << "," << before_y << ") "
	   << "-> (" << x_position.at(n) 
	   << "," <<  y_position.at(n) << ")."
	   << endl; **/
      
      //clean up a bit
      if(scale!=1){
	delete temp_single_int;
	delete weight_single_int;
      }


    }

    //if we are using the minimum error algorithm
    if(type==MINIMUM_ERROR){
      
      Complex_2D temp(lnx*scale,lny*scale);
      
      for(int i=0; i<lnx*scale; i++){
	for(int j=0; j<lny*scale; j++){
	  
	  int i_ = i - before_x*scale -  x_min*scale;
	  int j_ = j - before_y*scale -  y_min*scale;
	  
	  if( i_>=0 && j_>=0 && i_<nx && j_<ny ){
	    
	    temp.set_real(i,j,object->get_real(i_,j_));
	    temp.set_imag(i,j,object->get_imag(i_,j_));
	  }
	  
	}
      }
      
      object = &temp;
      
      nx = lnx*scale;
      ny = lny*scale;
      
      x_min = -x_position.at(n);
      y_min = -y_position.at(n);
      
      check_position(n,step_size, min_x, max_x, min_y, max_y);
      
       //return to normal
      object = &temp_object;
      nx = nx_c;
      ny = ny_c;
      x_min = x_min_c; 
      y_min = y_min_c;
    }

    cout << "moving prob " << n << " from ("
	 << before_x << "," << before_y << ") "
	 << "-> (" << x_position.at(n) 
	 << "," <<  y_position.at(n) << ")."
	 << endl;
   

    weights_set=false;
    add_to_object(n);
	
    if(forward)
      n++;
    else
      n--;
    
  }
  
  //set all the parameters back to normal
  object = pointer_to_object;

  nx = nx_c;
  ny = ny_c;
  x_min = x_min_c;
  y_min = y_min_c;
  parallel = parallel_c;
  beta=beta_c;

  //force the weights to be recalculated
  weights_set=false;


  //reset the object global transmission function
  //  scale_object(1-beta); 
  
  /**  for(int i=0; i<singleCDI.size(); i++){
    add_to_object(i);
    }**/


}

int PhaseDiverseCDI::check_position(int n_probe, double step_size, 
				    int min_x, int max_x,
				    int min_y, int max_y,
				    int tries){
  
  double x = x_position.at(n_probe);
  double y = y_position.at(n_probe);

  /**  cout << "checking probe "<< n_probe << " position, " 
       << x << "," << y << " with step size " 
       << step_size << ". Try no.: "<<tries<< endl;**/

  //done, we found the best position.
  if(step_size<1.0/scale)
    return SUCCESS;

  //failed, we moved around a bit, but couldn't find a local minima
  //in the error metric.

  if(tries > 10 || x < (x+min_x) || x > (x+max_x) || y < (y+min_y) || y > (y+max_y) ){
    cout << "Giving up on probe "<< n_probe << ". Could not find " 
	 << "it's position. Returning to the original. " 
	 << endl;
      return FAILURE;
  }
  
  BaseCDI * single = singleCDI.at(n_probe);
  
  double best_x=x;
  double best_y=y;
  double best_error=100;
   
  double new_x;
  double new_y;
  double size_x = single_result.at(n_probe)->get_size_x();
  double size_y = single_result.at(n_probe)->get_size_y();
  
  //try the 9 positions around the current one
  //record the one with the lowest error metric.

  //copy some local stuff.
  Complex_2D * object_copy = new Complex_2D(object->get_size_x(),object->get_size_y());
  object_copy->copy(*object);
  Complex_2D * single_copy = new Complex_2D(size_x,size_y);
  single_copy->copy(*single_result.at(n_probe));
  double beta_c = beta;

  beta = 0.5;
  weights_set=false;

  for(int i=-1; i<2; i++){
    for(int j=-1; j<2; j++){

      new_x=x+i*step_size;
      new_y=y+j*step_size;

      //checking whether the corrds are still within the 
      //boundary of the image
      //      if(new_x>0 && new_x<size_x && new_y>0 && new_y<size_y){

      x_position.at(n_probe)=new_x;
      y_position.at(n_probe)=new_y;

      add_to_object(n_probe);
      update_from_object(n_probe);

      /**      if(n_probe==4&&new_x==3294&&new_y==2917){
	static int h = 0;
	char blah[80];
	Double_2D temp(1024,1024);
	single_result.at(n_probe)->get_2d(PHASE,temp);
	sprintf(blah,"single_prob%i_try%i.tiff",n_probe,h);
	write_image(blah,temp,true,-0.1,0);

	object->get_2d(PHASE,temp);
	sprintf(blah,"object_prob%i_try%i.tiff",n_probe,h);
	write_image(blah,temp,true,-0.1,0);

	h++;

	cout << "error at ("<<new_x<<","<<new_y
	     << ")"<<" is "<<single->get_error()<<endl;

	     }**/

      single->iterate();
	
      if(single->get_error()<best_error){
	best_error = single->get_error();
	best_x = new_x;
	best_y = new_y;
      }

      //reset the object estimate:
      object->copy(*object_copy);
      //and
      set_result(single,*single_copy);
      single_result.at(n_probe)->copy(*single_copy);

    }
  }

  delete object_copy;
  delete single_copy;

  //set x and y to the best one :
  x_position.at(n_probe)=best_x;
  y_position.at(n_probe)=best_y;

  //recursively call this function with a smaller step size.
  if(best_x==x && best_y==y)
    step_size=step_size/2.0;
  
  int status = check_position(n_probe, step_size, 
			      min_x, max_x,
			      min_y, max_y,
			      ++tries);

  if(status == FAILURE){
    //return to orginal coordinates
    x_position.at(n_probe) = x ;
    y_position.at(n_probe) = y;
   return FAILURE;
  }

  /**  cout << "moving probe "<< n_probe << " by " 
       << x_position.at(n_probe)-x <<" in x "
       << "and "<< y_position.at(n_probe)-y<<" in y." 
       << endl; **/
  

  beta = 1.0;
  weights_set=false;

  return SUCCESS;
  
}


void PhaseDiverseCDI::set_up_weights(){
  
  if(weights_set)
    return;
  else
    weights_set=true;
  
  int frames = singleCDI.size();

  for(int n_probe=0; n_probe<frames; n_probe++){
    
    int lnx = single_result.at(n_probe)->get_size_x();
    int lny = single_result.at(n_probe)->get_size_y();
    BaseCDI * this_CDI = singleCDI.at(n_probe);
    double value;
    Double_2D illum_mag(lnx,lny);
    double max;

    if(typeid(*this_CDI)==typeid(FresnelCDI)){
      const Complex_2D & illum = ((FresnelCDI*)(this_CDI))->get_illumination_at_sample();
      illum.get_2d(MAG,illum_mag);
      max = illum_mag.get_max();
    }
    
    //    Double_2D support(lnx,lny);
    //this_CDI->get_support(support);
    
    for(int i_=0; i_< lnx; i_++){
      for(int j_=0; j_< lny; j_++){
	
	if(typeid(*this_CDI)==typeid(FresnelCDI))
	  value = illum_mag.get(i_,j_)/max;
	else
	  value = 1;

	value=alpha.at(n_probe)*pow(value,gamma);

	if(this_CDI->get_support().get(i_,j_)<=0.0)
	  value = 0;
	//else
	//  cout << support.get(i_,j_) <<endl;

	if(!parallel)
	  value*=beta;
	
	weights.at(n_probe)->set(i_,j_,value);

      }
    } 
  }
  
  if(parallel){
    
    Double_2D weight_norm(object->get_size_x(),object->get_size_y());
    
    for(int n=0; n<frames; n++){
      
      double x = x_position.at(n);
      double y = y_position.at(n);
      
      int this_size_x = single_result.at(n)->get_size_x();
      int this_size_y = single_result.at(n)->get_size_y();
      
      
      for(int i_=0; i_< this_size_x; i_++){
	for(int j_=0; j_< this_size_y; j_++){
	  
	  int i = get_global_x_pos(i_,x); // (i_-x-x_min)*scale;
	  int j = get_global_y_pos(j_,y); // (j_-y-y_min)*scale;
	  
	  if(i>=0&&j>=0&&i<nx&&j<ny){
	    
	    double new_weight = weights.at(n)->get(i_,j_);
	    double old_weight = weight_norm.get(i,j);
	    
	    if(n==0)
	      weight_norm.set(i,j,new_weight);
	    else
	      weight_norm.set(i,j,new_weight+old_weight);
	  }
	}
      }
    }
    
    for(int n=0; n<frames; n++){
      
      double x = x_position.at(n);
      double y = y_position.at(n);
      
      int this_size_x = single_result.at(n)->get_size_x();
      int this_size_y = single_result.at(n)->get_size_y();
      
      for(int i_=0; i_< this_size_x; i_++){
	for(int j_=0; j_< this_size_y; j_++){
	  
	  int i = get_global_x_pos(i_,x); 
	  int j = get_global_y_pos(j_,y); 

	  //	  int i = (i_-x-x_min)*scale;
	  // int j = (j_-y-y_min)*scale;
	  
	  if(i>=0&&j>=0&&i<nx&&j<ny){
	    
	    double old_weight = weights.at(n)->get(i_,j_);
	    double norm = weight_norm.get(i,j);

	    if(norm<=0)
	      weights.at(n)->set(i_,j_,0);
	    else
	      weights.at(n)->set(i_,j_,beta*old_weight/norm);

	  }
	}
      }
    } 
    //    write_image("weight_norm.tiff",weight_norm);
 
  }

  //  write_image("weight_0.tiff",*(weights.at(0)));
 };


  
/**
 * Update the result of the 'n_probe'th sub-frame to the
 * global object. This is a rather complicated, messy looking
 * function. But fast code has been prefered to clean code in this
 * case. 
 **/
void PhaseDiverseCDI::add_to_object(int n_probe){
  
  //get the result from the local FresnelCDI or PlanarCDI
  get_result(singleCDI.at(n_probe),*(single_result.at(n_probe)));  
  set_up_weights();
  
  double x_offset = x_position.at(n_probe);
  double y_offset = y_position.at(n_probe);
  
  //'small' holds the local frame result
  Complex_2D * small = single_result.at(n_probe);

  //'this_weight' holds the local frame weighting function
  Double_2D & this_weight = *(weights.at(n_probe));
  
  int lnx = small->get_size_x();
  int lny = small->get_size_y();
  
  //pretabulate the spacings (only used if sub-pixel reconstruction
  //is performed i.e. scale>1)
  double * sub_pos = new double[scale];
  for(int di=0; di < scale; di++){
    sub_pos[di] = (di + 0.5*((scale+1) % 2)) /((double) scale);
  }

  //i_, j_ - the local (small) coordinate system
  for(int i_=1; i_< lnx-1; i_++){
    for(int j_=1; j_< lny-1; j_++){
      
      double weight = this_weight.get(i_,j_);

      //don't even bother to do anything is the weight at this pixel is zero.
      if(weight!=0){

	//get the pixel positions in the global (object) frame.
	int i = get_global_x_pos(i_,x_offset)*scale;
	int j = get_global_y_pos(j_,y_offset)*scale;

	//bounds check
	if(i>=0&&j>=0&&i<nx&&j<ny){

	  //double sum = 0;

	  //get the value at the local pixel.
	  double f00r=small->get_real(i_,j_);
	  double f00i=small->get_imag(i_,j_);
	  
	  //if we are doing sub-pixel positioning, 
	  //(ie scale > 1) then we need to interpolate
	  //the 'small' array so it uses the same scale 
	  //as the object array.
	  if(scale!=1){

	    //get the values in the surrounding pixels
	    double f10r=small->get_real(i_+1,j_);
	    double f01r=small->get_real(i_,j_+1);
	    double f11r=small->get_real(i_+1,j_+1);
	    
	    double f10i=small->get_imag(i_+1,j_);
	    double f01i=small->get_imag(i_,j_+1);
	    double f11i=small->get_imag(i_+1,j_+1);
	    
	    //loop over the sub-pixel elements
	    for(int di=0; di < scale; di++){
	      for(int dj=0; dj < scale; dj++){

		double x = sub_pos[di];
		double y = sub_pos[dj];
		
		double x_1 = 1-x;
		double y_1 = 1-y;
		
		//work out the interpolated value for each
		//sub-pixel point
		double value_r = f00r*(x_1)*(y_1) + f10r*x*(y_1) 
		  + f01r*(x_1)*y + f11r*x*y;
		double value_i = f00i*(x_1)*(y_1) + f10i*x*(y_1) 
		  + f01i*(x_1)*y + f11i*x*y;
		
		//work out the global position for each sub-pixel point
		double new_i = (i_+0.5-x_offset-x_min)*scale + di;
		double new_j = (j_+0.5-y_offset-y_min)*scale + dj;
		
		//calculate the new value in the object (also using the pixel weight).
		double new_real = weight*value_r 
		  + object->get_real(new_i, new_j);
		double new_imag = weight*value_i 
		  + object->get_imag(new_i, new_j);	  
		
		//if in series mode...
		if(!parallel){
		  new_real -= weight*object->get_real(new_i, new_j);
		  new_imag -= weight*object->get_imag(new_i, new_j);
		}
		
		//set the final values
		object->set_real(new_i, new_j,new_real);
		object->set_imag(new_i, new_j,new_imag); 
	      
	      }
	    }
	  }
	  else{ //if we are not doing sub-pixel positioning

	    //just work out the value in the simple way.
	    double new_real = weight*f00r + object->get_real(i,j);
	    double new_imag = weight*f00i + object->get_imag(i,j);	  
	    
	    if(!parallel){
	      new_real -= weight*object->get_real(i, j);
	      new_imag -= weight*object->get_imag(i, j);
	    }
	    
	    object->set_real(i, j, new_real);
	    object->set_imag(i, j, new_imag); 

	  }

	}
      }
      
    }
  }

  /**  Double_2D temp(nx,ny);
  object->get_2d(MAG,temp);
  write_image("before_norm.tiff",temp,false,0,1);

  // normalise 
  if(scale!=1){
    for(int i_=1; i_< lnx-1; i_++){
      for(int j_=1; j_< lny-1; j_++){
      
	double weight = this_weight.get(i_,j_);
	
	//don't even bother to do anything is the weight at this pixel is zero.
	if(weight!=0){
	  
	  //get the pixel positions in the global (object) frame.
	  int i = get_global_x_pos(i_,x_offset)*scale;
	  int j = get_global_y_pos(j_,y_offset)*scale;
	  
	  //bounds check
	  if(i>=0&&j>=0&&i<nx&&j<ny){ 

	    double weight_r = 0;
	    double weight_i = 0;

	    for(int di=0; di < scale; di++){
	      for(int dj=0; dj < scale; dj++){
		//add to the sum for each pixel 
		//so we can normalise later
		weight_r+=object->get_real(i+di,j+dj);
		weight_i+=object->get_imag(i+di,j+dj);	
	      }
	    }
	
	    //	    cout << weight_i << " " <<  small->get_imag(i_,j_) << " ";
	    if(weight_r!=0&&weight_i!=0){
	      weight_r=small->get_real(i_,j_)*scale*scale/weight_r;
	      weight_i=small->get_imag(i_,j_)*scale*scale/weight_i;
	      
	      //	      cout << weight_i <<endl;
	      
	      for(int di=0; di < scale; di++){
		for(int dj=0; dj < scale; dj++){
		  object->set_real(i+di,j+dj,weight_r*object->get_real(i+di,j+dj));
		  object->set_imag(i+di,j+dj,weight_i*object->get_imag(i+di,j+dj));
		}
	      }   
	    }
	  }
	}
      }
    }
    } **/

  //  object->get_2d(MAG,temp);
  //write_image("after_norm.tiff",temp,false,0,1);

  delete [] sub_pos;
  
}

/**void PhaseDiverseCDI::get_object_sub_grid(Complex_2D & result,
					  double x_offset,
					  double y_offset){
  
  
    //using the local coordinate system
    for(int i_=0; i_< result.get_size_x(); i_++){
      for(int j_=0; j_< result.get_size_y(); j_++){
	
	//local coors (i,j) are global (object) coors  
	int i = i_ - x_offset*scale - x_min*scale;
	int j = j_ - y_offset*scale - y_min*scale;
	
	if(i>=0 && j>=0 && i<nx && j<ny){
	  result.set_real(i_, j_, object->get_real(i,j));
	  result.set_imag(i_, j_, object->get_imag(i,j));
	}
	else{
	  result.set_real(i_, j_, 1);
	  result.set_imag(i_, j_, 0);
	}
      }
    }
  
    };**/




//update the result from the 'n_probe'th frame using the
//global object function.
void PhaseDiverseCDI::update_from_object(int n_probe){

  set_up_weights();

  double x_offset = x_position.at(n_probe);
  double y_offset = y_position.at(n_probe);
  
  //small will hold a sub-array of the global object.
  Complex_2D * small = single_result.at(n_probe);

  Double_2D & this_weight = *weights.at(n_probe);

  int lnx = small->get_size_x();
  int lny = small->get_size_y();

  //i_,j_ - local coordinate system
  //i , j - global coordinate system
  for(int i_=0; i_< lnx; i_++){
      for(int j_=0; j_< lny; j_++){
	
	//only do something if the weight is non-zero. ie.
	//the 'n_probe'th frame contains information about this pixel.
	if(this_weight.get(i_,j_)!=0){

	  int i = get_global_x_pos(i_,x_offset)*scale;
	  int j = get_global_y_pos(j_,y_offset)*scale;

	  //bounds check
	  if(i>=0 && j>=0 && i<nx && j<ny){
	    
	    double real_value = 0;
	    double imag_value = 0;

	    //if we're doing sub-pixel reconstruction
	    //shrink the global (object) array down
	    //to the same dimensions as the local (single) one. 
	    if(scale!=1){
	      for(int di=0; di < scale; di++){
		for(int dj=0; dj < scale; dj++){
		  int new_i = i+di;
		  int new_j = j+dj;
		  real_value+=object->get_real(new_i,new_j);
		  imag_value+=object->get_imag(new_i,new_j);
		}
	      }
	      
	      double norm = 1/((double)(scale*scale));
	      //normalise to set the average (instead of the sum)
	      small->set_real(i_,j_,real_value*norm);
	      small->set_imag(i_,j_,imag_value*norm);
	    }
	    else{ //if we're not doing sub-pixel positioning, the
	          //just copy the global (object).
	      small->set_real(i_,j_,object->get_real(i,j));
	      small->set_imag(i_,j_,object->get_imag(i,j)); 
	    }
	  }
	}
      }
  }	    

  //set the result in the FresnelCDI or PlanarCDI object.
  set_result(singleCDI.at(n_probe),*(single_result.at(n_probe)));
  
};
 
