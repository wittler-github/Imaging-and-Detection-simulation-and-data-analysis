/*
 *  ELEtTTRA SImulation.c
 *  
 *
 *  Created by Michael Jones on 10/03/11.
 *  Copyright 2011 La Trobe University. All rights reserved.
 *
 */

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "io.h"
#include "Complex_2D.h"
#include "Double_2D.h"
#include "FresnelCDI.h"
#include <cstdlib>
#include <math.h>
#include <string>
#include <PlanarCDI.h>
//#include "TransmissionConstraint.h"

using namespace std;


/**********************************/

int main(int argc, char * argv[]){

//number of output iterations
  int output_iterations =50 ;
  int total_hio_iterations =200;
  int total_er_iterations = 200;
  int noise_level = 10;
  int shrink = 20; //apply shrinkwrap at this interval

  //array size
   int nx = 2048;
   int ny = 2048;
  //middle of the array
   int x_mid = nx/2;
   int y_mid = ny/2;

	//the field that passes through the slices, starts as white field, ends as ESW
	//not used for this simple example
	//Complex_2D wavefield(nx,ny); 

// experimental parameters and geometry
/********************************************/
	
       // double object_diameter = 8e-6;
	double lambda =4.891e-10;
	double zp_diameter = 160e-6; //zone plate diameter
	double zp_radius = zp_diameter/2.0;
	double zp_outer_zone = 50e-9;
	double beamstop = 20e-6; //radius of beamstop
	double focal = (zp_diameter * zp_outer_zone)/lambda; //zone plate focal length
	
	double zfs = 1.45e-3;// 18.513e-3 - 16.353e-3;// 0.6e-3;// 1.35e-3; //focal to sample
     	double zfd = 0.9;// 0.909513 -16.353e-3;// 0.5006;// 0.8; //focal to detector @ 0.8m, each pixel in sample plane is equal to: 18.75nm. 
	//for 2um object @ 18.75nm pixels, 106 pixels in object
	//for 3um object @ 18.75nm pixels, 160 pixels in object
 	//for 4um object @ 18.75nm pixels, 213 pixels in object
	//set relevant image stack folder: line 271
	double zsd =zfd-zfs;	
	double pix = 13.5e-6; //pixel size
	double zpr_f = zp_radius/focal;
	double far_field_radius = ((zp_radius/focal)*zfd)/pix; //beam radius in the detector
	int seed = 7;
	//copy the small into the large at your starting position
	int start_x = 640;
	int start_y = 350;
	double beamstop_radius = ((beamstop/focal)*zfs); //beamstop at sample
	double beam_stop_detector = ((beamstop/focal)*zfd)/pix; //beamstop at detector
	double beam_size_sample =  ((zp_radius/focal)*zfs);
	
 	double magnification = (zfs+zsd)/zfs;
        double object_diameter = (beam_size_sample - beamstop_radius);// ((nx/5)*pix)/magnification; 	
//run a few checks to see if the geometry is suitable
/****************************************************/
//will the sample fit in the beam?
      cout<<"beam size @ sample = "<<beam_size_sample<<endl;
	cout<<"beamstop @ sample = "<<beamstop_radius<<endl;
 	cout<<"beamstop @ detector(in pixels) ="<<beam_stop_detector<<endl;      
//zone plate diameter and outer zone
	cout<<"zone plate: diameter = "<<zp_diameter<<" , outer zone = "<<zp_outer_zone<<endl; 
//focal length
	cout<<"focal = "<<focal<<endl; 
 //distances
 cout<<"sample distance = "<<zfs<<" detector distance = "<<zfd<<endl;
 //will it fit on the detector
 cout<<"far_field_radius (in pixels) = "<<far_field_radius<<" needs to be less than "<<ny/2<<endl;
  	
  //is it sampled on the detector
 double sampled = (2.0 * object_diameter * pix)/lambda;
 cout<<"Sample to detector distance: "<<zsd<<" Must be larger than:"<<sampled<<endl;
 //what is the resolution
 double resolution = (lambda * zsd)/(nx * pix);
 cout<<"theoretical resolution = "<<resolution<<endl;
 //is the object sampled sufficiently
 double NF_sample = (object_diameter/2)*(object_diameter/2)/(lambda * zfs);
 cout<<"Fresnel at sample = "<<NF_sample<<endl;

 //check magnification and object diameter
cout << "Magnification =" <<magnification <<endl;
cout << "Object_diameter ="<< object_diameter << endl;

 /********************************************************/
	//create a circle to use as the whitefield support
		Double_2D circle(nx,ny); //array for the whitefield to be created
	double radius = far_field_radius;//use the radius of the beam defined by the geometry
	//create the circle using radius
	for(int y=y_mid-radius; y<=y_mid+radius; y++){
		for(int x=x_mid-radius; x<=x_mid+radius; x++){
			int y2 = abs(y - y_mid);
			int x2 = abs(x - x_mid);
			if(x2*x2+y2*y2 <= radius*radius)
				circle.set(x,y,1.0);
		}
	}
	
	//change from type double to int
	int beam_stop =beam_stop_detector; //size of the central stop
	
	//create the beamstop circle to subtract
	for(int y=y_mid-beam_stop; y<=y_mid+beam_stop; y++){
		for(int x=x_mid-beam_stop; x<=x_mid+beam_stop; x++){
			int y2 = abs(y - y_mid);
			int x2 = abs(x - x_mid);
			if(x2*x2+y2*y2 <= beam_stop*beam_stop)
				circle.set(x,y,0);
		}
	}
	
	//convolve it
	Complex_2D object_estimate2(nx,ny);
	PlanarCDI planar(object_estimate2,4);
	planar.convolve(circle,25,25); //25 & 25 gives the best results here

//convert Double_2D to complex_2D array for multiplying later
	Complex_2D circle_complex(nx,ny); 
	   for(int i=0; i<nx; i++){
		   for(int j=0; j<ny; j++){
			       circle_complex.set_value(i,j,REAL, 1.0/sqrt(2.0)*circle.get(i,j));
			       circle_complex.set_value(i,j,IMAG, 1.0/sqrt(2.0)*circle.get(i,j));
			     }
		   }
					
//create the phase and magnitude of incident whitefield at the detector 
	Complex_2D incident_whitefield_det(nx,ny); 
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			
			//distance (r) from the centre of the array to each individual pixel
		double rho_2_d =sqrt( pow(pix*(i-x_mid),2.0) + pow(pix*(j-y_mid),2.0));
		
		//phase delay at the detector
		//double phi_B_d = -M_PI*( (pow(rho_2_d,2) )/(lambda*zfd) );		
		double phi_B_d =((M_PI*(pow(rho_2_d,2.0)))/(lambda))*((-1.0/((zfd-zfs)))+(1.0/(zfd)));
		//amplitude
		int size = 600; //size defines the size of the gaussian spot
		double amp = exp ( -( rho_2_d)/(2 * size * pix));
			
		incident_whitefield_det.set_real(i,j,amp * cos(phi_B_d));
		incident_whitefield_det.set_imag(i,j,amp * sin(phi_B_d));
		
				
	}
}
	
	//multiply by the size of the beam in the detector 
	// as defined by the circle and geometry
	incident_whitefield_det.multiply(circle_complex);
	
	//set outside white field to zero
	for( int i =0; i<nx; i++){
		for( int j = 0; j<ny; j++){
			if(circle.get(i,j)==0){//set outside to zero
				incident_whitefield_det.set_real(i,j,0);
				incident_whitefield_det.set_imag(i,j,0);
			}
		}
	}
	
//write the whitefield to tiffs to check
Double_2D result_wf1(nx,ny);
Double_2D result_wf2(nx,ny);
ostringstream phase_check ( ostringstream::out );
phase_check << "Fresnel_NAS_pics/phase_check.tiff";
ostringstream mag_check ( ostringstream::out );
mag_check << "Fresnel_NAS_pics/magnitude_check.tiff";
incident_whitefield_det.get_2d(PHASE,result_wf1);
write_tiff(phase_check.str(),result_wf1,true);
incident_whitefield_det.get_2d(MAG,result_wf2);
write_tiff(mag_check.str(),result_wf2);

		
	/***********************************
   This projection will propagate the whitefield at the detector to the sample plane
  *********************************************************/
  //make a 2D array and allocate some memory.
  //This will be used to output the image of the 
  //current estimate.

  Double_2D result(nx,ny);
  //copy the whitefield to wf. this way incident_whitefield_det is still at the detector plane, and wf will be at the sample plane
  Complex_2D wf(nx,ny);
  wf.copy(incident_whitefield_det);  


  //set parameters
  double wavelength = lambda; //wavelength
  double fd_wf = zfd; //focal to detector
  double fs_wf = zfs; //focal to sample
  double ps = pix; //pixel size

  //create the projection object which will be used later to
  //perform the reconstuction.
  Complex_2D object_estimate(nx,ny);

  FresnelCDI whitefield(object_estimate, //estimate
		  wf, //white field 
		  wavelength, //wavelength
		  fd_wf, //focal-to-detector
		  fs_wf, //focal-to-sample
		  ps); //pixel size

  //******************************************
  //Project the object to the image plane 
  //******************************************  
  
   //get the illuminating white field in the sample plane
  whitefield.propagate_from_detector(wf);

  //check it
  Double_2D wf_sample(nx,ny);
  wf.get_2d(MAG_SQ,wf_sample);
  write_tiff("Fresnel_NAS_pics/wf_sample.tiff",wf_sample);
 //check it
  Double_2D phase_sample(nx,ny);
  wf.get_2d(PHASE,phase_sample);
  write_tiff("Fresnel_NAS_pics/phase_sample.tiff",phase_sample);

								/***********************use a pointer here rather than writing the .cplx********/
  //write to a .cplx file
  write_cplx("Fresnel_NAS_pics/whitefield_at_start.cplx",wf);
  write_cplx("Fresnel_NAS_pics/wf_at_detector.cplx",incident_whitefield_det);	
	
}
