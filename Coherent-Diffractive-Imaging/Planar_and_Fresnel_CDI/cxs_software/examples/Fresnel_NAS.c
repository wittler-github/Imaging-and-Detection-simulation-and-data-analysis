/**
 * @file Fresnel_NAS.c 
*/
 

#include <iostream>
#include <math.h>
#include <string>
#include <cstdlib> 
#include "io.h"
#include "Complex_2D.h"
#include "Double_2D.h"
//#include "PlanarCDI.h"
#include "FresnelCDI.h"
#include <sstream>

#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <typeinfo>

using namespace std;


/**************************************/
int main(void){

  //File to write error metric, SNR in
  ofstream output_SNR_EM_C("Fresnel_NAS_pics/SNR_EM/SNR_EM_C/SNR_EM.txt");
  ofstream output_SNR_EM_M("Fresnel_NAS_pics/SNR_EM/SNR_EM_M/SNR_EM.txt");
  ofstream output_SNR_EM_P("Fresnel_NAS_pics/SNR_EM/SNR_EM_P/SNR_EM.txt");
  ofstream output("Fresnel_NAS_pics/F_T_1.txt");//TEST
  ofstream output2("Fresnel_NAS_pics/F_T_2.txt");//TEST
  ofstream output5("Fresnel_NAS_pics/F_T_5.txt");//TEST
  ofstream output6("Fresnel_NAS_pics/F_T_6.txt");//TEST
  ofstream Rec_and_True_TF_MAX_MIN("Fresnel_NAS_pics/Rec_&_True_TF_MAX_MIN.txt");
  

  //Define some quantities which will be used in the code:

  //the object file name
  const static char * object_file_name = "Fresnel_NAS_pics/CXS_Article_image_2048_s_GB_15pxhv.tif";//Fresnel_NAS_pics/FCDI_HEI_2048_8.tiff";//FCDI_HEI_2048_8_smooth.tiff"CXS_2048.tif";//

  //the file which provides the support (pixels with the value 0
  //are considered as outside the object)
  const static char * support_file_name = "Fresnel_NAS_pics/FCDI_HEI_2048_8_wf_SUPPORT.tiff";//Fresnel_NAS_pics/FCDI_HEI_2048_8_Rectangle_Support.tiff";//FCDI_HEI_2048_8_wf_SUPPORT.tiff";//FCDI_HEI_2048_8_SHAPE_SUPPORT.tiff";//;//FCDI_HEI_2048_16_support.tiff";

  //number of hybrid input-out iterations to perform.
  const int hio_iterations = 52;

  //the number of error-reduction iterations to perform.
  const int er_iterations = 52;

  //output the current image ever "output_iterations"
  const int output_iterations = 13;


  /****** get the object from an image file ****************/

  //get the object from file
  
  Double_2D object;

  //read the data into an array
  int status = read_tiff(object_file_name, object);
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }

  int n_x = object.get_size_x();
  int n_y = object.get_size_y();


  /*****Get the support from file**********/
  Double_2D support(n_x,n_y);
  status = read_tiff(support_file_name, support);

      //make sure the object and support are the same dimensions
  if(support.get_size_x() != n_x || support.get_size_y() != n_y){
    cout << "dimensions do not match....exiting" <<endl;
    return(1);
  }

  /***Get the simulated wf in detector plane****/
    Complex_2D wfdata(n_x,n_y);
    status = read_cplx("Fresnel_NAS_pics/wf_at_detector.cplx", wfdata);
  if(!status){ //give an error if we couldn't open the file.
    cout  << "Maybe you need to run ./FCDI_WF_example.exe "
	  << " or WF_Fresnel_NAS first... exiting"  << endl;
    return(1);
  }


//--------Reconstruct phase for the whitefield----//
Double_2D WF_T_R(n_x,n_y);
wfdata.get_2d(MAG_SQ,WF_T_R);
 
  status = write_dbin("Fresnel_NAS_pics/WF_TR_MAG_SQ.dbin", WF_T_R);
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }

  status = system("./FresnelCDI_WF_example.exe");


  status = read_cplx("Fresnel_NAS_pics/wf_recovered.cplx", wfdata);
  if(!status){ //give an error if we couldn't open the file.
    cout  << "Maybe you need to run ./FCDI_WF_example.exe "
	  << ".. exiting"  << endl;
    return(1);
  }
//-----------------------------------------------------------------------------------------//


  //For use to write results to images
  Double_2D result(n_x,n_y);

  /**** Create the projection/reconstruction object *****/

  //set the experimental parameters (all are in meters)
  double wavelength = 4.891e-10;// 4.892e-10; //wavelength
  double fd = 0.9; //focal to detector
  double fs = 1.45e-3; //focal to sample
  double ps = 13.5e-6; //pixel size

  //create the projection object which will be used later to
  //perform the reconstuction.
  Complex_2D first_guess(n_x,n_y);

  FresnelCDI my_fresnel(first_guess, //estimate
		  wfdata, //white field 
		  wavelength, //wavelength
		  fd, //focal-to-detector
		  fs, //focal-to-sample
		  ps, //pixel size
		  1.0); //normalisation of white-field to sample data

  //***************************************************
  //  Simulate a diffraction pattern
  //***************************************************  

  /*** First we will create a transmission function **/

  // create the complex array which it will be stored in 
  Complex_2D input(n_x,n_y);

  //lets say our sample is 150nm maximum thickness and is made from gold
  double thickness = 150e-9;
  double max = object.get_max();
  cout << "MAX" << max << endl;
  
  int th_C = 0.2*max;
  int th_Mg = 0.5*max;
  int th_Au = 0.3*max;
  
  

  double d_C = 7.246e-5;
  double b_C = 1.283e-6;

  double d_Mg = 5.711e-5;
  double b_Mg = 7.274e-6;

  double d_Au = 3.211e-4; //2.8713e-4; // 6.45e-4; 
  double b_Au = 1.888e-4; //2.2334e-4; // 1.43e-4;

 



  double k = (2.0 * M_PI)/wavelength;

  double mag = 0;
  double phase = 0;

  double pix_th = 0;

  for (int i = 0; i<n_x; i++){
    for (int j = 0; j<n_y; j++){

      //calculate the amplitude and phase of the transmission function
      //at each point (from the image).
      pix_th = thickness*object.get(i,j)/max;
      
    
      if(i >= 690 && i <= 910 && j >= 650 && j <= 990){
        mag = exp(-k*b_C*pix_th);
	phase = -k*d_C*pix_th;
	input.set_value(i,j,MAG,mag);
	input.set_value(i,j,PHASE,phase);
        }
	
      if(i >= 910 && i <= 1135 && j >= 650 && j <= 990){
        mag = exp(-k*b_Mg*pix_th);
	phase = -k*d_Mg*pix_th;
	input.set_value(i,j,MAG,mag);
	input.set_value(i,j,PHASE,phase);
        }
	
      if(i >= 1135 && i <= 1360 && j >= 650 && j <= 990){
        mag = exp(-k*b_Au*pix_th);
	phase = -k*d_Au*pix_th;
	input.set_value(i,j,MAG,mag);
	input.set_value(i,j,PHASE,phase);
        }

        mag = exp(-k*b_C*pix_th);
	phase = -k*d_C*pix_th;
	input.set_value(i,j,MAG,mag);
	input.set_value(i,j,PHASE,phase);

    }
  }
      

  write_tiff("Fresnel_NAS_pics/OBJECT_C_X_S.tiff",object);
  cout << "OBJECT_C_X_S.tiff" << endl;

  //Save the magnitude and phase to file so we know how it looks
  input.get_2d(MAG,result);
  Rec_and_True_TF_MAX_MIN << "MAX true_t_f mag " << result.get_max() << endl;
  Rec_and_True_TF_MAX_MIN << "MIN true_t_f mag " << result.get_min() << endl;
  write_tiff("Fresnel_NAS_pics/trans_mag_input.tiff",result);

  input.get_2d(PHASE,result);
  Rec_and_True_TF_MAX_MIN << "MAX true_t_f phase " << result.get_max() << endl;
  Rec_and_True_TF_MAX_MIN << "MIN true_t_f phase " << result.get_min() << endl;
  write_tiff("Fresnel_NAS_pics/trans_phase_input.tiff",result);

  Complex_2D true_t_f(n_x,n_y);
  true_t_f.copy(input);
 



  /********* Now simulate: ***********************************/ 
   
  //get the illuminating white field in the sample plane
  my_fresnel.propagate_from_detector(wfdata);

  //TEST write wfdata//
  Double_2D TESTWFDATA(n_x,n_y);
  wfdata.get_2d(MAG,TESTWFDATA);
  write_tiff("Fresnel_NAS_pics/TESTWFDATA",TESTWFDATA);
  //END TEST//

  //multiply the transmission function by the white field
  input.multiply(wfdata);

  //get the magnitude of the wave
  input.get_2d(MAG,result);
  //write the output to file
  write_tiff("Fresnel_NAS_pics/sample_plane.tiff",result);


  //get the transmission function and write it
  Complex_2D tfc(n_x,n_y);
  Double_2D tfd_MAG(n_x,n_y);
  Double_2D tfd_PH(n_x,n_y);
  my_fresnel.get_the_transmission_function(input,tfc,true);
  tfc.get_2d(MAG, tfd_MAG);
  tfc.get_2d(PHASE, tfd_PH);
  write_tiff("Fresnel_NAS_pics/tfd_MAG.tiff",tfd_MAG);
  write_tiff("Fresnel_NAS_pics/tfd_PH.tiff",tfd_PH);
  

  //propagate to the detector
  my_fresnel.propagate_to_detector(input);

  //take the wf back to the detector so we're ready for the reconstructions
  my_fresnel.propagate_to_detector(wfdata);

  //write input intensity to Double_2D array
  Double_2D intensitydata(n_x,n_y); 
  input.get_2d(MAG_SQ,intensitydata);

  //write the diffraction pattern to a file
  write_tiff("Fresnel_NAS_pics/forward_projection.tiff",intensitydata,true);
  
  //change relaxation parameter
  my_fresnel.set_relaxation_parameter(1);  

//-----DEFINED QUANTITIES USED IN FOR LOOPS, IN ORDER----------//

  //Used for photon scaling, gain and SNR..
  double nop_ps = 1e+7;//white-field flux photons per second
  double nop_array[11] = {0.5e-1*nop_ps,1e-1*nop_ps,5e-1*nop_ps,1e+1*nop_ps,1e+2*nop_ps,5e+2*nop_ps,1e+3*nop_ps,5e+3*nop_ps,1e+4*nop_ps,1e+5*nop_ps,1e+7*nop_ps};  
  double nop = 0;//number of photons present in the measurement
  double si = 0;//sum all elements in the array
  double aipp = 0; //average intensity per photon
  double gain = 157;//gain of ccd ADU's per photon  
  double SNR = 0; //to be used in S/N calculation
  double sum_Ik=0,sum_srIk=0,sum_Bk=0; 

  //Used in poisson noise
  double mu = 0, wf_mu = 0;
  unsigned int p = 0, wf_p = 0;

  //Used in converting nr of photons to ADU's
  double m = 0, wf_m = 0;  
  double sd_ADC = 30;
 
  //Used in bias addition
  double ps_array[9] = {0,0.1,0.3,0.5,0.7,0.9,0.95,0.99,1};
  double a_bl = 0; 
  double wf_a_bl = 0;
  double n_o_f = 100;//nr of 1s frames
  double meanb = 210;//mean bias level for gaussian distribution
  double sdb = 14;// ;//sd bias level in ADU for gaussian distribution
  Double_2D mbb(n_x,n_y);
  double wf_in_int = 0;

  //Array to reset intensity & whitefield for each loop
  Double_2D intensity(n_x,n_y);
  Complex_2D wf(n_x,n_y);

  //Used for true object
  Complex_2D True(n_x,n_y);
  Double_2D true_image(n_x,n_y);

  //To be used when sorting the 5 best reconstructions//
  Complex_2D ER_5_b_array[5] = { Complex_2D(n_x,n_y), Complex_2D(n_x,n_y), 
                   Complex_2D(n_x,n_y), Complex_2D(n_x,n_y), Complex_2D(n_x,n_y) };
  
  Complex_2D HIO_5_b_array[5] = { Complex_2D(n_x,n_y), Complex_2D(n_x,n_y), 
                   Complex_2D(n_x,n_y), Complex_2D(n_x,n_y), Complex_2D(n_x,n_y) };

  //used when writing out the 5 best reconstructions
  Double_2D TEST(n_x,n_y);

  //Used to write white-field
  Complex_2D WF_R(n_x,n_y);
  Complex_2D TESTWFC(n_x,n_y);
  Double_2D TESTWFD(n_x,n_y),WFTEST(n_x,n_y);
  Double_2D WFD(n_x,n_y);
  Double_2D WF_TR(n_x,n_y);

  //Used for qualatative tests of code 
  Double_2D TESTER(n_x,n_y);
  Double_2D TESTHIO(n_x,n_y);

  //Used for writing the reconstructions det.amp and phase
  Complex_2D ER1DPC(n_x,n_y);
  Complex_2D HIO1DPC(n_x,n_y);
  Double_2D ER_1_AMP_det_ls(n_x,n_y), ER_1_PH_sam(n_x,n_y);
  Double_2D HIO_1_AMP_det_ls(n_x,n_y), HIO_1_PH_sam(n_x,n_y);

  //To be used when writing reconstructed quantities to tiff
  Complex_2D cu(n_x,n_y);
  Complex_2D TFu(n_x,n_y);
  Double_2D m_or_p(n_x,n_y);

//-----FOR LOOP DIFFERENT PHOTON NUMBER------------------//
for( int pn = 10; pn <=10  ; pn++){
  
//------FOR LOOP FOR BIAS NOISE--//
//for( int bl = 0 ; bl<= 8 ; bl++ ){

  intensity.copy(intensitydata);
  
  wf.copy(wfdata);

  cout << "HEJSVEJS" << endl;    
  nop = nop_array[pn];//number of photons present in the measurement
  si = intensity.get_sum();//sum all elements in the array
  aipp = si/nop; //average intensity per photon
  
  //Reset SNR
  SNR = 0; 
  sum_Ik=0,sum_srIk=0,sum_Bk=0; 
  

  //----To scale the simulated diffraction data along with wf----// 
  intensity.scale( 1/aipp );
  wf.scale(sqrt(1/aipp));
  
   
//----------Generate true image-------------------------//
  True.copy(input);//set equal to complex diffraction pattern at detector
  
  True.scale( sqrt(gain/aipp) );//convert intensity to corrsponding number of photons
    
  True.add(wf,-sqrt(gain));
  
  my_fresnel.propagate_from_detector(True);//propagate wavefield to sample plane 
  True.get_2d(MAG,true_image);
  write_tiff("Fresnel_NAS_pics/Trueimage.tiff",true_image,false);

  
  //---Check the true transmission function------
  Complex_2D TTT(n_x,n_y);
  Double_2D TTTM(n_x,n_y),TTTP(n_x,n_y);
  my_fresnel.get_the_transmission_function(True,TTT,true);
  TTT.get_2d(MAG,TTTM);
  TTT.get_2d(PHASE,TTTP);
  write_tiff("Fresnel_NAS_pics/TTTM",TTTM,false);
  write_tiff("Fresnel_NAS_pics/TTTP",TTTP,false);
  //-------------------------------------------
  


  //Set up for random number generation--------------------------
  gsl_rng *r1,*r2,*r3,*r4,*r5,*r6;

  r1 = gsl_rng_alloc (gsl_rng_mt19937);
  r2 = gsl_rng_alloc (gsl_rng_mt19937);
  r3 = gsl_rng_alloc (gsl_rng_mt19937);
  r4 = gsl_rng_alloc (gsl_rng_mt19937);
  r5 = gsl_rng_alloc (gsl_rng_mt19937);
  r6 = gsl_rng_alloc (gsl_rng_mt19937);

  gsl_rng_set ( r1, 15);
  gsl_rng_set ( r2, 25);
  gsl_rng_set ( r3, 35);
  gsl_rng_set ( r4, 45);
  gsl_rng_set ( r5, 55);
  gsl_rng_set ( r6, 65);
  //-------------------------------------------------------------


 //------Apply inherent poisson random noise-----------------//

  for(int i=0;i <n_x;i++){
    for(int j=0;j <n_y;j++){
	 
	mu = round( intensity.get(i,j) ); //mean of poisson distribution	
	
	p = gsl_ran_poisson (r1, mu); //random variate from poisson distribution with mean mu	
	     
	intensity.set(i,j,p);
 
	
        wf_mu = round( wf.get_value(i,j,MAG_SQ) );
 
        wf_p = gsl_ran_poisson (r2, wf_mu );
 
        wf.set_value(i,j,MAG,sqrt(wf_p) );
	
        
	
     }	
   }

 //--Convert nr of photons to ADU'S--	
  
  for(int i=0;i< n_x;i++){
    for(int j=0;j< n_y;j++){
	  
	m = intensity.get(i,j)*gain + gsl_ran_gaussian(r3,sqrt(intensity.get(i,j)*pow(sd_ADC,2)) );//Equals a sum of gaussian random variates with respective mean=gain and sd 30

	if(m < 0){
	 m = 0;
	}

	m = round(m);
	         
        sum_Ik += m;
        sum_srIk += sqrt(m);

	intensity.set(i,j,m);
	
	
        //-Below same for whitefield--//
	
        wf_m = wf.get_value(i,j,MAG_SQ)*gain + gsl_ran_gaussian(r4,sqrt(wf.get_value(i,j,MAG_SQ)*pow(sd_ADC,2)) );
 
        if(wf_m < 0){
 	  wf_m = 0;
 	}

        wf_m = round(wf_m);
		
	wf.set_value(i,j,MAG,sqrt(wf_m));
                          
    }
  }


/*
//-------Gaussian distributed bias level with already subtracted constant for the background-------/
	
  for(int i=0;i<n_x;i++){
    for(int j=0;j<n_y;j++){
	
        a_bl = n_o_f*meanb + gsl_ran_gaussian(r5, sqrt(n_o_f*pow(sdb,2)) );
               
        a_bl -= n_o_f*meanb*ps_array[bl];
        
	if( a_bl < 0){
	   a_bl = 0;
	}
        
	sum_Bk += a_bl;
		
	mbb.set(i,j,a_bl);//Matrice to add to reciprocal space

	//Same for white-field
        wf_a_bl = n_o_f*meanb + gsl_ran_gaussian(r6, sqrt(n_o_f*pow(sdb,2)) );
       
 	wf_a_bl -= n_o_f*meanb*ps_array[bl];

	if( wf_a_bl < 0){
	   wf_a_bl = 0;
	}

	wf_in_int = wf.get_value(i,j,MAG_SQ);//initial intensity

	wf.set_value(i,j,MAG, sqrt(wf_in_int + wf_a_bl));//initial intensity + added bias

      }
   }

  intensity.add(mbb,1);


//---------------------------------------------------------------
*/

// free random number generators
  gsl_rng_free (r1);
  gsl_rng_free (r2);
  gsl_rng_free (r3);
  gsl_rng_free (r4);
  gsl_rng_free (r5);
  gsl_rng_free (r6);
 

  SNR = sum_Ik/(sum_srIk + sum_Bk);
  printf(" SNR %f",SNR);


  //write the output to file
  write_tiff("Fresnel_NAS_pics/sim_intensity.tiff",intensity,false);
  write_tiff("Fresnel_NAS_pics/sim_intensity_ls.tiff",intensity,true);

 
  //Reset normalisation of white-field to sample dp at detector
   Double_2D WF_intensity(n_x,n_y);
   wf.get_2d(MAG_SQ,WF_intensity);
   double norm_sdp_to_wf = ( intensity.get_sum() )/( WF_intensity.get_sum() );
   cout << "NORM SDP TO WF" << norm_sdp_to_wf<< endl;
   //my_fresnel.set_normalisation(norm_sdp_to_wf);

   //--------Reconstruct phase for the noise added whitefield----//
Double_2D WF_N_T_R(n_x,n_y);
wf.get_2d(MAG_SQ,WF_N_T_R);
 
  status = write_dbin("Fresnel_NAS_pics/WF_TR_MAG_SQ.dbin", WF_N_T_R);
  if(!status){
    cout << "failed.. exiting"  << endl;
    return(1);
  }

  status = system("./FresnelCDI_WF_example.exe");


  status = read_cplx("Fresnel_NAS_pics/wf_recovered.cplx", wf);
  if(!status){ //give an error if we couldn't open the file.
    cout  << "Maybe you need to run ./FCDI_WF_example.exe "
	  << ".. exiting"  << endl;
    return(1);
  }
//-----------------------------------------------------------------------------------------//
 

  //--Reset the whitefield to be used in reconstruction--//
   my_fresnel.set_illumination(wf);

   /*
   //--TEST
   my_fresnel.get_illumination(TESTWFC);
   TESTWFC.get_2d(MAG_SQ,TESTWFD);
   wf.get_2d(MAG_SQ,WFTEST);
   write_tiff("Fresnel_NAS_pics/TEST_ILL", TESTWFD,false );
   write_tiff("Fresnel_NAS_pics/TEST_ILL_ls", TESTWFD,true );
   write_tiff("Fresnel_NAS_pics/TEST_WF", WFTEST,false );
   write_tiff("Fresnel_NAS_pics/TEST_WF_ls", WFTEST,true );
   //--TEST
   */
  

  

  
 //Reset quantities used in sorting the 5 best rec.
  double ER_5_b_e_array [5] = { 5,5,5,5,5 };
  double HIO_5_b_e_array [5] = { 5,5,5,5,5 };
  int place = 0;


//***FOR LOOP 5 RANDOM STARTS***//
for( int l = 0; l < 5 ; l++ ){
 

  /*************** Do the reconstruction *******************/
  
  //-------FOR HIO------------------//
if( hio_iterations > 0 ){
  
  my_fresnel.set_support(support);
  
  my_fresnel.set_intensity(intensity);
  
  my_fresnel.set_algorithm(HIO);

  //set the initial guess to be random inside the support
  //and zero outside. Note that this must be called
  //after "my_fresnel.set_support()"
  my_fresnel.initialise_estimate(6*l);
  
  //apply the projection operators
  for(int i=0; i <= hio_iterations; i++){

    cout << "iteration " << i << endl;

    my_fresnel.iterate();

    if(i%output_iterations==0){

      ostringstream temp_str ( ostringstream::out ) ;
      first_guess.get_2d(MAG,result);
      temp_str << "Fresnel_NAS_pics/HIO_sim_result_" << i << ".ppm";
      write_ppm(temp_str.str(), result);

      //my_fresnel.apply_shrinkwrap(1.5,0.12);
     
    }
  }
 
  //---SORT INTO ARRAYS---//
  place = 0;  
  for(  ; place < 5 && my_fresnel.get_error() > HIO_5_b_e_array[place]; place++); 
 
  //we found a new best estimate
  if(5>0 && place < 5){

    Complex_2D * temp = new Complex_2D(n_x,n_y);
  
    for(int i=4; i>place; i--){
      HIO_5_b_e_array[i] = HIO_5_b_e_array[i-1];
      HIO_5_b_array[i].copy( HIO_5_b_array[i-1] );
    }
    
    HIO_5_b_e_array[place] = my_fresnel.get_error();
  
    temp->copy(first_guess);
    
    HIO_5_b_array[place].copy(*temp);      

    delete temp;
  } 
 

}
//-------FOR HIO-END-------------//

//-------FOR ER----------------//
if( er_iterations > 0 ){
  
  my_fresnel.set_support(support);

  my_fresnel.set_intensity(intensity);  

  my_fresnel.set_algorithm(ER);

  my_fresnel.initialise_estimate(6*l);

  for(int i=0; i <= er_iterations; i++){
    
    cout << "iteration " << i << endl;
    
    my_fresnel.iterate();

    if(i%output_iterations==0){

      ostringstream temp_str ( ostringstream::out ) ;
      first_guess.get_2d(MAG,result);
      temp_str << "Fresnel_NAS_pics/ER_sim_result_" << i << ".ppm";
      write_ppm(temp_str.str(), result);

      //my_fresnel.apply_shrinkwrap(1.5,0.12);
      
    }
  }


 //-----SORT INTO ARRAYS-----//
  place = 0;  
  for(  ; place < 5 && my_fresnel.get_error() > ER_5_b_e_array[place]; place++); 
 
  //we found a new best estimate
  if(5>0 && place < 5){

    Complex_2D * temp = new Complex_2D(n_x,n_y);
  
    for(int i=4; i>place; i--){
      ER_5_b_e_array[i] = ER_5_b_e_array[i-1];
      ER_5_b_array[i].copy( ER_5_b_array[i-1] );
    }
    
    ER_5_b_e_array[place] = my_fresnel.get_error();
  
    temp->copy(first_guess);
    
    ER_5_b_array[place].copy(*temp);      

    delete temp;
  } 
 
 }

//-------FOR ER--END-------------------//
  

  printf( "ERROR: %f \n", my_fresnel.get_error() );//TEST
  printf(" SNR %f \n",SNR);//TEST

}//for loop random ESW

//----WRITE TIFFS OF PHOTON NUMBER RESULTS-----//


//-Write out the 5 best reconstruction---//

for( int i =0; i<5 ; i++){

   ostringstream temp_str ( ostringstream::out ) ;
   HIO_5_b_array[i].get_2d(MAG,TEST);
   temp_str << "Fresnel_NAS_pics/PN/PN_" << pn << "/HIO_5_" << i << ".ppm";
   write_ppm(temp_str.str(), TEST);
   output << HIO_5_b_e_array[i] << endl;

}

for( int i =0; i<5 ; i++){

   ostringstream temp_str ( ostringstream::out ) ;
   ER_5_b_array[i].get_2d(MAG,TEST);
   temp_str << "Fresnel_NAS_pics/PN/PN_" << pn << "/ER_5_" << i << ".ppm";
   write_ppm(temp_str.str(), TEST);
   output << ER_5_b_e_array[i] << endl;

}
//---------------//



//Write tiff's of simulated wf and intensity
   ostringstream t_s_Noisy_DP_AMP_ls( ostringstream::out );
   t_s_Noisy_DP_AMP_ls << "Fresnel_NAS_pics/PN/PN_"<< pn <<"/Noisy_DP_AMP_ls.tiff";
   write_tiff(t_s_Noisy_DP_AMP_ls.str(), intensity,true);
  
   ostringstream t_s_Noisy_WF_AMP_ls( ostringstream::out );
   wf.get_2d(MAG_SQ,WFD);
   t_s_Noisy_WF_AMP_ls << "Fresnel_NAS_pics/PN/PN_"<< pn <<"/Noisy_WF_AMP_ls.tiff";
   write_tiff(t_s_Noisy_WF_AMP_ls.str(), WFD,true);


//Write diffraction pattern and phase for initial and 1'st ER, HIO
   ER1DPC.copy( ER_5_b_array[0] );
   HIO1DPC.copy( HIO_5_b_array[0] );

   ER1DPC.get_2d(PHASE, ER_1_PH_sam);
   ostringstream t_s_ER_1_PH_sam( ostringstream::out );
   t_s_ER_1_PH_sam << "Fresnel_NAS_pics/PN/PN_"<< pn <<"/ER_1_PH_sam.tiff";
   write_tiff(t_s_ER_1_PH_sam.str(), ER_1_PH_sam);

   HIO1DPC.get_2d(PHASE, HIO_1_PH_sam);
   ostringstream t_s_HIO_1_PH_sam( ostringstream::out );
   t_s_HIO_1_PH_sam << "Fresnel_NAS_pics/PN/PN_"<< pn <<"/HIO_1_PH_sam.tiff";
   write_tiff(t_s_HIO_1_PH_sam.str(), HIO_1_PH_sam);

   my_fresnel.propagate_to_detector( ER1DPC );
   my_fresnel.propagate_to_detector( HIO1DPC );

   ER1DPC.get_2d(MAG_SQ, ER_1_AMP_det_ls);
   ostringstream t_s_ER_1_AMP_det_ls( ostringstream::out );
   t_s_ER_1_AMP_det_ls << "Fresnel_NAS_pics/PN/PN_"<< pn <<"/ER_1_AMP_det_ls.tiff";
   write_tiff(t_s_ER_1_AMP_det_ls.str(), ER_1_AMP_det_ls, true);

   HIO1DPC.get_2d(MAG_SQ, HIO_1_AMP_det_ls);
   ostringstream t_s_HIO_1_AMP_det_ls( ostringstream::out );
   t_s_HIO_1_AMP_det_ls << "Fresnel_NAS_pics/PN/PN_"<< pn <<"/HIO_1_AMP_det_ls.tiff";
   write_tiff(t_s_HIO_1_AMP_det_ls.str(), HIO_1_AMP_det_ls, true);


//Output the max and min of the ER1, HIO1 phase to know how the tiff scales from -pi to +pi
   output5 << ER_1_PH_sam.get_min() << setw(15) << ER_1_PH_sam.get_max() << setw(15) << HIO_1_PH_sam.get_min() << setw(15) << HIO_1_PH_sam.get_max() << endl;



//-----Write reconstructions----//

//nr 1 mod ER mag
cu.copy(ER_5_b_array[0]);
my_fresnel.app_modulus(cu);
cu.get_2d(MAG,m_or_p);
ostringstream t_s_ER_1_MAG_sam_mod( ostringstream::out );
t_s_ER_1_MAG_sam_mod << "Fresnel_NAS_pics/PN/PN_"<< pn <<"/ER_1_MAG_sam_mod.tiff";
write_tiff(t_s_ER_1_MAG_sam_mod.str(), m_or_p);
//nr 1 mod ER phase
cu.copy(ER_5_b_array[0]);
my_fresnel.app_modulus(cu);
cu.get_2d(PHASE,m_or_p);
ostringstream t_s_ER_1_PH_sam_mod( ostringstream::out );
t_s_ER_1_PH_sam_mod << "Fresnel_NAS_pics/PN/PN_"<< pn <<"/ER_1_PH_sam_mod.tiff";
write_tiff(t_s_ER_1_PH_sam_mod.str(), m_or_p);

//nr 1 mod TF ER mag
my_fresnel.get_the_transmission_function(cu, TFu, 1);
TFu.get_2d(MAG,m_or_p);
Rec_and_True_TF_MAX_MIN << "MAX 1 mod TF ER mag " << m_or_p.get_max() << endl;
Rec_and_True_TF_MAX_MIN << "MIN nr 1 mod TF ER mag " << m_or_p.get_min() << endl;
ostringstream t_s_TF_ER_1_MAG_sam_mod( ostringstream::out );
t_s_TF_ER_1_MAG_sam_mod << "Fresnel_NAS_pics/PN/PN_"<< pn <<"/TF_ER_1_MAG_sam_mod.tiff";
write_tiff(t_s_TF_ER_1_MAG_sam_mod.str(), m_or_p);
//nr 1 mod TF ER phase
TFu.get_2d(PHASE,m_or_p);
Rec_and_True_TF_MAX_MIN << "MAX 1 mod TF ER phase " << m_or_p.get_max() << endl;
Rec_and_True_TF_MAX_MIN << "MIN nr 1 mod TF ER phase " << m_or_p.get_min() << endl;
ostringstream t_s_TF_ER_1_PH_sam_mod( ostringstream::out );
t_s_TF_ER_1_PH_sam_mod << "Fresnel_NAS_pics/PN/PN_"<< pn <<"/TF_ER_1_PH_sam_mod.tiff";
write_tiff(t_s_TF_ER_1_PH_sam_mod.str(), m_or_p);

//nr 1 mod HIO mag
cu.copy(HIO_5_b_array[0]);
my_fresnel.app_modulus(cu);
cu.get_2d(MAG,m_or_p);
ostringstream t_s_HIO_1_MAG_sam_mod( ostringstream::out );
t_s_HIO_1_MAG_sam_mod << "Fresnel_NAS_pics/PN/PN_"<< pn <<"/HIO_1_MAG_sam_mod.tiff";
write_tiff(t_s_HIO_1_MAG_sam_mod.str(), m_or_p);
//nr 1 mod HIO phase
cu.copy(HIO_5_b_array[0]);
my_fresnel.app_modulus(cu);
cu.get_2d(PHASE,m_or_p);
ostringstream t_s_HIO_1_PH_sam_mod( ostringstream::out );
t_s_HIO_1_PH_sam_mod << "Fresnel_NAS_pics/PN/PN_"<< pn <<"/HIO_1_PH_sam_mod.tiff";
write_tiff(t_s_HIO_1_PH_sam_mod.str(), m_or_p);

//nr 1 mod TF HIO mag
my_fresnel.get_the_transmission_function(cu, TFu, 1);
TFu.get_2d(MAG,m_or_p);
ostringstream t_s_TF_HIO_1_MAG_sam_mod( ostringstream::out );
t_s_TF_HIO_1_MAG_sam_mod << "Fresnel_NAS_pics/PN/PN_"<< pn <<"/TF_HIO_1_MAG_sam_mod.tiff";
write_tiff(t_s_TF_HIO_1_MAG_sam_mod.str(), m_or_p);
//nr 1 mod TF HIO phase
TFu.get_2d(PHASE,m_or_p);
ostringstream t_s_TF_HIO_1_PH_sam_mod( ostringstream::out );
t_s_TF_HIO_1_PH_sam_mod << "Fresnel_NAS_pics/PN/PN_"<< pn <<"/TF_HIO_1_PH_sam_mod.tiff";
write_tiff(t_s_TF_HIO_1_PH_sam_mod.str(), m_or_p);

//-----END WRITE RESULTS PHOTON NUMBER-----//


/////////////////////////////////////////////////////////
//----COMPLEX ESW EM-----//


double em10ER = my_fresnel.C_f_em( ER_5_b_array[0], true, True, false);

double em12ER = my_fresnel.C_f_em( ER_5_b_array[0], true, ER_5_b_array[1], true);

double csemER = my_fresnel.C_f_em( ER_5_b_array[0], false, ER_5_b_array[0], true);


double em10HIO = my_fresnel.C_f_em( HIO_5_b_array[0], true, True, false);

double em12HIO = my_fresnel.C_f_em( HIO_5_b_array[0], true, HIO_5_b_array[1], true);

double csemHIO = my_fresnel.C_f_em( HIO_5_b_array[0], false, HIO_5_b_array[0], true);


//....COMPLEX TF EM...//
double T_em10ER = my_fresnel.C_O_F_T_EM( ER_5_b_array[0], true, true_t_f, 0, 0);

double T_em12ER = my_fresnel.C_O_F_T_EM( ER_5_b_array[0], true, ER_5_b_array[1], true, 1);

double T_csemER = my_fresnel.C_O_F_T_EM( ER_5_b_array[0], false, ER_5_b_array[0], true, 1);


double T_em10HIO = my_fresnel.C_O_F_T_EM( HIO_5_b_array[0], true, true_t_f, 0, 0);

double T_em12HIO = my_fresnel.C_O_F_T_EM( HIO_5_b_array[0], true, HIO_5_b_array[1], true, 1);

double T_csemHIO = my_fresnel.C_O_F_T_EM( HIO_5_b_array[0], false, HIO_5_b_array[0], true, 1);


output_SNR_EM_C << setw(15) << SNR << setw(15) << em10ER << setw(15) << em12ER << setw(15) << csemER << setw(15)
                                      << em10HIO << setw(15) << em12HIO << setw(15) << csemHIO << setw(15) << T_em10ER << setw(15) << T_em12ER << setw(15) << T_csemER << setw(15)
                                      << T_em10HIO << setw(15) << T_em12HIO << setw(15) << T_csemHIO <<  endl;
////////////////////////////////////////////////////////////


}//for loop photon nr or bias 

     
  

  return 0;
}

