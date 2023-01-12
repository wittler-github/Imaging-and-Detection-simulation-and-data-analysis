/**
 * @file PhaseDiverse_example.c
 *
 * This example file shows how phase-diverse / ptychographic
 * reconstruction can be performed. The example uses data from Corey
 * Putkunz which can be found on osiris at
 * /data/cputkunz/phase_diverse_cdi/example_data.tar.gz
 * 
 *
 * @author Nadia Davidson <nadiamd@unimelb.edu.au> 
 *
 * Last modified on 15/4/2011
 *
 */

//this one for the GAMS for APS AUG11

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "io.h"
#include "TransmissionConstraint.h"
#include "Complex_2D.h"
#include "Double_2D.h"
#include "FresnelCDI.h"
#include "PhaseDiverseCDI.h"
#include <cstdlib>
#include <math.h>
#include <string>
#include "utils.h"

using namespace std;

//transmission constraint, phase <=0, |T| <=1;
void my_charge_flipping(Complex_2D & transmission){
 
  double min_phase = -M_PI/2;
  double max_phase = 0;
  double min_T = 0;
  double max_T =1;

 //we know the Transmission can't be out by more than
  for(int i=0; i<transmission.get_size_x(); i++){
    for(int j=0; j<transmission.get_size_y(); j++){

      if(transmission.get_mag(i,j)>max_T)
	transmission.set_mag(i,j,max_T);

      if(transmission.get_mag(i,j)<min_T)
	transmission.set_mag(i,j,min_T);
    }
  }
  //we know the phase can't be out by more than
  for(int i=0; i<transmission.get_size_x(); i++){
    for(int j=0; j<transmission.get_size_y(); j++){

      if(transmission.get_phase(i,j)>max_phase)
	transmission.set_phase(i,j,max_phase);

      if(transmission.get_phase(i,j)<min_phase)
	transmission.set_phase(i,j,min_phase);
    }
  }

}

/**********************************/

int main(int argc, char * argv[]){

  ////////////////////////
  int nx = 2048;
  int ny = 2048;

  const int max_iterations = 100;
  const int out_image = 50;

  //make the support
  Double_2D beam(nx,ny);
  double beam_fraction = 0.5;
  double i0 = (nx-1)/2.0;
  double j0 = (ny-1)/2.0;
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      if(sqrt((i-i0)*(i-i0)+(j-j0)*(j-j0)) 
	 < beam_fraction*sqrt(j0*j0+i0*i0))
	beam.set(i,j,100);
      else{
	beam.set(i,j,0);
      }
    }
  }

  //Make a new PhaseDiverseCDI object
  PhaseDiverseCDI pd(1.0,
		     3.0,
		     true,
		     1);

  //how many frames make up the reconstruction
  const int frames = 1;

  //initialise some arrays to hold data.
  FresnelCDI * proj[frames];
  Complex_2D * object_estimate[frames];
  Complex_2D * wf;
  Double_2D * diffraction;
  Double_2D * this_support;

  //set the experimental parameters (all are in meters)
  double wavelength = 4.892e-10; //wavelength
  double fz = 9.814e-3; //zone plate to focal distance
  double ps = 13.5e-6; //pixel size

 //zone plate to detector distances
  double fd = 0.83+fz;
		    
  //sample to zone plate distances
  double zfs = 1.0550e-3;
  double zfs2 = zfs - 100e-6;
  double fs[frames] = {
		       fz+zfs, //I
		       fz+zfs, //J
		       fz+zfs, //M
		       fz+zfs, //N
		      
		       fz+zfs2, //I2
		       fz+zfs2, //J2
			   fz+zfs2, //M2
		       fz+zfs2 //N2
		       
};
  
  //white-field normalisation
  double norm[frames] = {
	  0.0221136,
	  0.0219296,
	  0.0220714,
	  0.0221987,
	  0.0221136,
	  0.0219296,
	  0.0220714,
	  0.0221987
  }; 
  


  
  //reconstructed white-field file names
  std::string wf_string[frames] = {				   
				
				"/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane1/WF_MEMBRANE/wf_recovered_mem.cplx",
				"/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane1/WF_MEMBRANE/wf_recovered_mem.cplx",
				"/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane1/WF_MEMBRANE/wf_recovered_mem.cplx",
				"/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane1/WF_MEMBRANE/wf_recovered_mem.cplx",
				"/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane2/WF_MEM2/wf_recovered_mem.cplx",
				"/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane2/WF_MEM2/wf_recovered_mem.cplx",
				"/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane2/WF_MEM2/wf_recovered_mem.cplx",
				"/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane2/WF_MEM2/wf_recovered_mem.cplx"}; 

  //data file names
  
  std::string diff_string[frames] = {
					"/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane1/POS_I/I_1secST.dbin",
				    "/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane1/POS_J/J_1secST.dbin",
   				    "/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane1/POS_M/M_1secST.dbin",
				    "/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane1/POS_N/N_1secST.dbin",
		
				    "/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane2/I2/I2_1secST.dbin",
				    "/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane2/J2/J2_1secST.dbin",
				    "/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane2/M2/M2_1secST.dbin",
				    "/media/DATA1/APS_AUG_2011/APS_AUG11_GAM_PTYCHO/GAMS_TIFF/Plane2/N2/N2_1secST.dbin"
					   
};
  
  
//these are the coordinates in sample plane pixels given a 1.5 um spacing
//use the formula;
// DeltaX = (lambda * ZD)/(NPixels * PixelSize)
//		  = (4.9e-10 * 0.83)/(2048 * 13.5e-6)
//		  =~ 15 nm
// So; 1.5 um spacing equals 1500 nm / 15 nm =~ 100 pixels. There are of course slight errors in reality as you can see below.
// In the theroetical model just make it 100 pixels.
	
  double x_pos[frames] = {
			  0,
			  100,
			  0,
			  100,
			  
			
			  -2,
			  103,
			  -2,
			  103,
			 
}; 
  
  double y_pos[frames] = {
			  141,
			  141,
			  232,
			  232,
			 
			  162,
			  162,
			  260,
			  260,
			  
}; 
  TransmissionConstraint tc4;
  // tc4.set_enforce_unity(true);
  //tc4.set_charge_flipping(false);
  tc4.set_custom_constraint(my_charge_flipping);
  
  for(int n=0; n < frames; n++){

    //Set up the fresnel CDI in the same way you would
    //if you weren't doing phase-diversity
    object_estimate[n] = new Complex_2D(nx,ny);

    //read the white-field from file
    wf = new Complex_2D(nx,ny);
    read_cplx(wf_string[n], *wf);
    
    //read the data from file
    diffraction = new Double_2D(nx,ny);
    read_dbin(diff_string[n],nx,ny,*diffraction);
    
    //set-up the reconstruction for a single frame
    proj[n]= new FresnelCDI(*object_estimate[n],
			    *wf,
			    wavelength,
			    fd-fz,
			    fs[n]-fz,
			    ps,
			    norm[n]);

    //set the support and intensity and initialise
    proj[n]->set_intensity(*diffraction);
    proj[n]->set_support(beam);
    proj[n]->initialise_estimate();

    //add a complex constraint
   proj[n]->set_complex_constraint(tc4);

    delete wf;
    delete diffraction;

    //New part.. Add the FresnelCDI to the PhaseDiverseCDI.
    pd.add_new_position(proj[n], x_pos[n], y_pos[n]);

  }

  //initialise the phase diverse transmission function
  pd.initialise_estimate();
  
  //now run the reconstruction
  for(int i=0; i< max_iterations; i++){

    //do some position alignment
    //if(i==0)
    //   pd.adjust_positions(PhaseDiverseCDI::CROSS_CORRELATION);
    //  if(i==10)
    //    pd.adjust_positions(PhaseDiverseCDI::MINIMUM_ERROR,true,-10,10,-10,10);
    
    //iteration!
    pd.iterate();
     if(i%out_image==0)
      {

     Complex_2D * object = pd.get_transmission();
ostringstream temp_str ( ostringstream::out ) ;
	  temp_str << "object_phase_GAMS_1expx1_" << i << ".cplx";
      write_cplx(temp_str.str(),*object);
      temp_str.clear();
      }
	   }


  Complex_2D * object = pd.get_transmission();

  Double_2D result(object->get_size_x(),object->get_size_y());

  //write the magnitude and phase of it to an image file
  // object->get_2d(MAG,result);
  write_image("object_mag_recovered_1expx1.tiff",result,false);

  object->get_2d(PHASE,result);
  write_image("object_phase_recovered_1expx1.tiff",result,false, -1.2,0);
  
  //save the transmission function in case we want to use it later.
  write_cplx("trans_1expx2.cplx",*object);


  return 0;
  
}
