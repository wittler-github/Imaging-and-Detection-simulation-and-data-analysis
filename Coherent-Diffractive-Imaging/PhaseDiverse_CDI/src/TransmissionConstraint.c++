// Copyright 2011 Nadia Davidson for The ARC Centre of Excellence in 
// Coherent X-ray Science. This program is distributed under the GNU  
// General Public License. We also ask that you cite this software in 
// publications where you made use of it for any part of the data     
// analysis. 

#include <cstdlib>
#include <iostream>
#include <Double_2D.h>
#include <Complex_2D.h>
#include <TransmissionConstraint.h>



/**
 * @param phase_sign The sign of the phase to be fliped. ie. If you
 * pass the function 1 is will flip all positive phases to
 * negative. If you give it -1, it will flip all negative phases to
 * positive.
 *
 */

using namespace std;
  
TransmissionConstraint::TransmissionConstraint(){
    region_map = NULL; 
    custom_constraint = NULL;

    set_enforce_unity(true);
    set_charge_flipping(true);

  }

TransmissionConstraint::~TransmissionConstraint(){
  if(region_map!=0){
    delete region_map;
  }
}
  
void TransmissionConstraint::add_complex_constraint(ComplexConstraint & new_constraint){

  //  cout << "Adding ComplexConstraint region"<<endl;

    int nx = new_constraint.get_region()->get_size_x();
    int ny = new_constraint.get_region()->get_size_y();

    if(region_map==0){
      region_map = new Double_2D(nx,ny);
    }

    if(region_map->get_size_x() != nx || region_map->get_size_y() != ny){
      cout << "Dimensions of the complex constraint region "
	   << "do not match the others given previously. "
	   << "Ignoring this complex constraint."<<endl;
    }
    
    complex_constraint_list.push_back(&new_constraint);

    int pos =  complex_constraint_list.size();
    
    for(int i=0; i < nx; i++){
      for(int j=0; j < ny; j++){
	if(new_constraint.get_region()->get(i,j)>0){
	  if(region_map->get(i,j)>0){
	    cout <<"Overlap between complex constraint region"
		 <<" detected. The last constraint to be set "
		 << "will be ignored in"
		 << "the overlapping region."<<endl;
	  }
	  else
	    region_map->set(i,j,pos);
	}
      }
    }
}


  /** This is the fundamental method in this class */
void TransmissionConstraint::apply_constraint(Complex_2D & transmission){

    int nx = transmission.get_size_x();
    int ny = transmission.get_size_y();

    int region_number;
    ComplexConstraint * current_constraint = 0;

    vector<double> mean_lnA;
    vector<double> mean_phase;
    
    int regions = complex_constraint_list.size();

    if(regions>0){

      //reset the mean c values to zero
      for(int i=0; i < regions; i++){
	mean_lnA.push_back(0);
	mean_phase.push_back(0);
      }

      //recalculate the mean c for each region
      for(int i=0; i < nx; i++){
	for(int j=0; j < ny; j++){
	
	  if(region_map->get(i,j)>0){
	    region_number = region_map->get(i,j)-1;
	    mean_phase.at(region_number) += fabs(transmission.get_phase(i,j));
	    mean_lnA.at(region_number) += fabs(log(transmission.get_mag(i,j)));
	  }
	}
      }

      for(int i=0; i<regions; i++){
	current_constraint = complex_constraint_list.at(i);
	current_constraint->set_c_mean( mean_lnA.at(i) / mean_phase.at(i) );
      }
    }
    
    /**  now we start the main loop to alter the
	   values in the transmission function array */
    
    double mag_old;
    double mag_new;
    double phase_old;
    double phase_new;

    if(regions>0||do_charge_flip||do_enforce_unity){

      for(int i=0; i < nx; i++){
	for(int j=0; j < ny; j++){
	
	  phase_old = transmission.get_phase(i,j);
	  mag_old = transmission.get_mag(i,j);

	  //apply complex constraints based on refractive indicies.
	  if(regions>0 && region_map->get(i,j)>0){
	    region_number =  region_map->get(i,j)-1;
	    current_constraint = complex_constraint_list.at(region_number);
	    
	    //cout <<"mag/phase old: "<< phase_old<<" "<< mag_old<<endl;

	    mag_new = current_constraint->get_new_mag(mag_old, phase_old);
	    phase_new = current_constraint->get_new_phase(mag_old, phase_old);
	    
	    //	cout <<"mag/phase old: "<< phase_new<<" "<< mag_new<<endl;

	    transmission.set_mag(i,j,mag_new);
	    transmission.set_phase(i,j,phase_new);
	     
	    mag_old=mag_new;
	    //phase_old=transmission.get_phase(i,j);

	  }
	  
	  if(do_charge_flip && flip_sign*transmission.get_imag(i,j)>0)
	    transmission.set_imag(i,j,-1*transmission.get_imag(i,j));
      

	  //restrict the transmission function to unit if requested
	  if(do_enforce_unity && mag_old>1)
	    transmission.set_mag(i,j,1);

	}
      }
    }

    //finally, run a custom constraint if one has been set!.
    if(custom_constraint)
      (*custom_constraint)(transmission);

}


