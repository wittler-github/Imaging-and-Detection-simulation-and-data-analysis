//author:  Nadia Davidson <nadiamd@unimelb.edu.au>
//date last modified: 12/1/2011

/**
 * @file CDI_reconstruction.c
 *
 * \a CDI_reconstruction.exe A tool for performing Planar or Fresnel
 * ESW reconstruction.  This tool is provided as a demonstrative tool
 * and to obtain results quickly.
 *
 * \par Usage: CDI_reconstruction.exe \<config filename\> \<reco_type\> \<seed\>
 *
 * 
 * where reco_type may be:
 * - "planar" - planar
 * - "fresnel_wf" - fresnel white-field reconstruction (3-plane propagation) 
 * - "fresnel" - fresnel object reconstruction. with the white-field 
 *               previously reconstructed.
 *
 * The seed should be an integer. If the seed is excluded from the
 * command line arguments, it is assumed to be "0". If reco_type is
 * also excluded, it is assumed to be "planar".
 *
 * \par Example:
 * \verbatim CDI_reconstruction.exe planar_example.config "planar" 3 \endverbatim
 * Perform planar CDI reconstruction using the configuration given in the file,
 * "planar_example.config". The random number generator (used to initial the
 * first guess) is given a seed value of 3.
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <cstdlib> 
#include "io.h"
#include "Complex_2D.h"
#include "Double_2D.h"
#include "PlanarCDI.h"
#include "FresnelCDI.h"
#include "FresnelCDI_WF.h"
#include "Config.h"

using namespace std;

#define OUTPUT_MINIMAL 0
#define OUTPUT_ITER    1
#define OUTPUT_ERROR   2

static const string planar_string="planar";
static const string fresnel_string="fresnel";
static const string fresnel_wf_string="fresnel_wf";


void print_usage(){

  cout << "Usage: " << endl << endl
       << "CDI_reconstruction.exe <config filename> " 
       << "<reco_type> <seed>" << endl << endl
       << "where <reco_type> may be: " << planar_string 
       << ", " << fresnel_string
       << " or " << fresnel_wf_string << endl
       << "<seed> should be an integer" << endl <<endl
       << "If <reco_type> and <seed> do not need to be specified."
       << "if they are not, <reco_type> will default to "<<planar_string 
       << " and <seed> to 0." << endl;

}

int main(int argc, char * argv[]){

  /** work out which config file to use **/
  string config_file = "";

  //and set the seed of the initial guess
  int seed = 0;

  string reco_type = "";

  if(argc==1){
    cout << endl << "No config file given, ";
    print_usage();
    exit(0);
  }
  else
    config_file=argv[1];

  if(argc<=2){
    cout << "No reconstruction type given... using planar" <<endl;
    reco_type = planar_string;
  }
  else
    reco_type = argv[2];

  if(argc<=3)
    cout << "No seed value given... using the default (0)" <<endl;
  else
    seed = atoi(argv[3]);

  if(argc>4){
    cout << endl << "Wrong number of arguments given. ";
    print_usage();
    exit(0);
  }

  Config c(config_file);
  if(c.getStatus()==FAILURE){
    cout << "Could not open the file "<< config_file<< ". Exiting" <<endl;
    exit(0);
  }

  /** read the config file **/

  //output the current image every "output_iterations"
  int output_iterations = c.getDouble("output_iterations");

  //names of the algorithms to use in order
  list<string> * algorithms = c.getStringList("algorithm");
 
  //the number of iterations to perform for each algorithm
  list<int> * iterations = c.getIntList("iterations");

  //get the number of pixels in x and y of the image
  const int pixels_x = c.getInt("pixels_x");
  const int pixels_y = c.getInt("pixels_y");
  
  //do some error checking. Were all the values we need present
  //in the config file?
  if(c.getStatus()==FAILURE){
    cout << "Problem reading the configuration file. Exiting" <<endl;
    exit(0);
  }
  if(algorithms->size()!=iterations->size()){
    cout << "The number of algorithms and the number "
	 <<"of iterations do not match. Exiting"<< endl;
    exit(0);
  }

  /** get some optional configuration parameters **/
  int output_level = c.getInt("info_level");
  int output_diffraction_estimate = c.getInt("output_diffraction_estimate");
  int use_log_scale_for_diffraction = c.getInt("use_log_scale_for_diffraction");
  int use_log_scale_for_object = c.getInt("use_log_scale_for_object");

  int shrinkwrap_iterations = c.getInt("shrinkwrap_iterations");
  double shrinkwrap_gauss_width = c.getDouble("shrinkwrap_gauss_width");
  double shrinkwrap_threshold = c.getDouble("shrinkwrap_threshold");

  string output_file_type = c.getString("output_file_type");

  /*******  set up the reconstruction ***************/

  //create the projection object which will be used to
  //perform the reconstruction.
  Complex_2D object_estimate(pixels_x,pixels_y);

  string starting_point_file_name = c.getString("starting_point_file_name");

  //if a file name has been given try to load it from the file    
  if(starting_point_file_name.compare("")!=0){
    if(!read_cplx(starting_point_file_name,object_estimate)){
      cout << "Can not process the file "<< starting_point_file_name 
	   << ".. exiting"  << endl;
      return(1);
    }
  }

  PlanarCDI * proj = 0;
  
  //the data file name
  string data_file_name = c.getString("data_file_name");

  //the file which provides the support (pixels with the value 0
  //are considered as outside the object)
  string support_file_name = c.getString("support_file_name");

  //output filename prefix
  string output_file_name_prefix = c.getString("output_file_name_prefix");
  ostringstream temp_str ( ostringstream::out ) ;
  temp_str << output_file_name_prefix << ".cplx" << flush;
  string output_file_name = temp_str.str();

  if(reco_type.compare(planar_string)==0){ //if Planar CDI
    proj = new PlanarCDI(object_estimate);
  }
  else{
    /** get the experimental parameters needed for FCDI */
    double beam_wavelength = c.getDouble("beam_wavelength");
    double zone_focal_length = c.getDouble("zone_focal_length");
    double focal_detector_length = c.getDouble("focal_detector_length");
    double focal_sample_length = c.getDouble("focal_sample_length");
    double pixel_size = c.getDouble("pixel_size");
    double normalisation = c.getDouble("normalisation");
    

    if(reco_type.compare(fresnel_string)==0){ //if Fresnel CDI

      //open the file where the reconstructed white field is
      Complex_2D white_field(pixels_x,pixels_y);
      int status = read_cplx(c.getString("white_field_reco_file_name"), 
			     white_field); 
      //check that the file could be opened okay
      if(!status){
	cout << "failed to read white-field data " 
	     <<".. exiting"  << endl;
	return(1);
      }

      cout << white_field.get_value(0,0,MAG_SQ) << endl;

      //create the iterator object
      proj = new FresnelCDI(object_estimate,
			    white_field, 
			    beam_wavelength, 
			    focal_detector_length, 
			    focal_sample_length, 
			    pixel_size, 
			    normalisation);
      
    }
    //if reconstructing the Fresnel white field.
    else if(reco_type.compare(fresnel_wf_string)==0){

      proj = new FresnelCDI_WF(object_estimate,
			       beam_wavelength, 
			       zone_focal_length, 
			       focal_detector_length, 
			       pixel_size);
      output_file_name = c.getString("white_field_reco_file_name");
      data_file_name = c.getString("white_field_data_file_name");
      support_file_name = c.getString("white_field_support_file_name");

      //reset the algorithm and number of iterations.
      algorithms->clear();
      algorithms->push_back("ER"); //this is a dummy value
      //since only one algorithm is available for white-field reco.
      iterations->clear();
      iterations->push_back(c.getInt("wf_iterations"));
      //don't use shrink-wrap either
      shrinkwrap_iterations = 0; 
    }
    else{
      cout << "the reconstruction type specified ("
	   << reco_type << ") is "
	   << "unrecognised. Please use: "
	   << planar_string << ", "
	   << fresnel_string << " or "
	   << fresnel_wf_string << endl;
      exit(1);
    }
  }

  /*** get the diffraction data from file and read into an array ***/
  Double_2D data;
  read_image(data_file_name, data, pixels_x, pixels_y);  
  
  //read_image does the error checking for us and would exit if the file
  //was not read

  if( pixels_x != data.get_size_x() || pixels_y != data.get_size_y() ){
    cout << "Dimensions of the data to not match those given ... exiting"  << endl;
    return(1);
  }

  /******* get the support from file and read it into an array *****/
  Double_2D support;
  read_image(support_file_name, support, pixels_x, pixels_y);  
  if( pixels_x != support.get_size_x() || pixels_y != support.get_size_y() ){
    cout << "Dimensions of the support to not match ... exiting"  << endl;
    return(1);
  }

  //set the support and intensity
  proj->set_support(support);
  proj->set_intensity(data);

  //Initialise the current object ESW with a random numbers
  if(starting_point_file_name.compare("")==0)
    proj->initialise_estimate(seed);

  //make a 2D array and allocate some memory.
  //This will be used to output the image of the 
  //current estimate.
  Double_2D result(pixels_x,pixels_y);

  /******* run the reconstruction *********/

  list<string>::iterator algorithms_itr = algorithms->begin();
  list<int>::iterator iterations_itr = iterations->begin();

  //loop over the algorithms
  int i=0;
  int cumulative_iterations = 0;
  while(algorithms_itr != algorithms->end()&&
	iterations_itr != iterations->end()){

    if(output_level!=OUTPUT_MINIMAL)
      cout << "Switching to the "<< (*algorithms_itr)
    	   <<" algorithm" << endl;
    
    //get the projection
    int alg = PlanarCDI::getAlgFromName(*algorithms_itr);
    if(alg == -1 ){
      std::cout << "Could not find reconstruction algorithm"
		<< " with the name "<< (*algorithms_itr)
		<< ". Exiting" << endl;
      exit(0);
    }

    proj->set_algorithm(alg);
    cumulative_iterations+=(*iterations_itr);
    
    for(; i < cumulative_iterations; i++){
      if(output_level!=OUTPUT_MINIMAL)
	cout << "Iteration " << i << endl;
      

      //apply the iterations  
      proj->iterate(); 
      if(output_level==OUTPUT_ERROR && i>0)
	cout << "Error for iteration "<<(i-1)<<" is " 
	     << proj->get_error() << endl;

      if(i%output_iterations==0){
	//output the current estimate of the object
	ostringstream temp_str ( ostringstream::out ) ;
	object_estimate.get_2d(MAG,result);
	temp_str << output_file_name_prefix << "_" << i << "."
		 << output_file_type << flush;
	write_image(temp_str.str(), result, use_log_scale_for_object);
	//temp_str.clear();

	//output the estimation of the intensity in 
	//the detector plane if needed
	if(output_diffraction_estimate){
	  Complex_2D * temp = object_estimate.clone();
	  proj->propagate_to_detector(*temp);
	  temp->get_2d(MAG_SQ,result);
	  temp_str << output_file_name_prefix 
		   << "_diffraction_" << i 
		   << "."<< output_file_type << flush;
	  write_image(temp_str.str(), result, 
		    use_log_scale_for_diffraction); 
	  delete temp;
	}
      }
      if(shrinkwrap_iterations!=0&&
	 i%shrinkwrap_iterations==(shrinkwrap_iterations-1)){
	proj->apply_shrinkwrap(shrinkwrap_gauss_width, 
			       shrinkwrap_threshold);
	
      }
    }
    iterations_itr++;
    algorithms_itr++;
  }

  //write out the final result
  write_cplx(output_file_name, object_estimate);

  //if it's fresnel reconstruction also output the
  //transmission function
  if(reco_type.compare(fresnel_string)==0){
    ((FresnelCDI*) proj)->get_transmission_function(object_estimate);
    write_cplx("transmission_function.cplx", object_estimate);
  }

  
  //clean up
  delete algorithms;
  delete iterations;
  delete proj;
  
  return 0;
}

