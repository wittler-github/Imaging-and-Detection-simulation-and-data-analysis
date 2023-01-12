#include <iostream>
#include <fstream>
//#include <stdlib.h>
#include <vector>
#include <sstream>
#include <cmath>

#if defined(HAVE_HDF_MFHDF_H)
#include "hdf/mfhdf.h"
#else
#include "mfhdf.h"
#endif

#include "io.h"
#include "Double_2D.h"

using namespace std;

#define FAILURE 0
#define SUCCESS 1

/**
 * @class hdf_anonymous_array
 * 
 * This is bad code which allows me to get around the fact that
 * I don't know the number of bytes per pixel in the data file 
 * until it's opened.
 */
class hdf_anonymous_array{

  int type_;
  int8 * a_int8;
  uint8  * a_uint8;
  int16  * a_int16;
  uint16 * a_uint16;
  int32 * a_int32;
  uint32 * a_uint32;
  float32 * a_float32;
  float64 * a_float64;
  
 public:

  hdf_anonymous_array(int type, int size){
    type_ = type;

    a_int8=0;
    a_uint8=0;
    a_int16=0;
    a_uint16=0;
    a_int32=0;
    a_uint32=0;
    a_float32=0;
    a_float64=0;

    switch( type ){
    case DFNT_INT8:
      a_int8 = new int8[size];
      break;
    case DFNT_UINT8:
      a_uint8 = new uint8[size];
      break;
    case DFNT_INT16:
      a_int16 = new int16[size];
      break;
    case DFNT_UINT16:
      a_uint16 = new uint16[size];
      break;
    case DFNT_INT32:
      a_int32 = new int32[size];
      break;
    case DFNT_UINT32:
      a_uint32 = new uint32[size];
      break;
    case DFNT_FLOAT32:
      a_float32 = new float32[size];
      break;
    case DFNT_FLOAT64:
      a_float64 = new float64[size];
      break;
      
    default:
      cout << "Not familiar with the HDF4 type.. exiting.." << endl;
    }    
  }
  
  ~hdf_anonymous_array(){
    switch( type_ ){
    case DFNT_INT8:
      delete[] a_int8;
      break;
    case DFNT_UINT8:
      delete[] a_uint8;
      break;
    case DFNT_INT16:
      delete[] a_int16;
      break;
    case DFNT_UINT16:
      delete[] a_uint16;
      break;
    case DFNT_INT32:
      delete[] a_int32;
      break;
    case DFNT_UINT32:
      delete[] a_uint32;
      break;
    case DFNT_FLOAT32:
      delete[] a_float32;
      break;
    case DFNT_FLOAT64:
      delete[] a_float64;
      break;
      
    default:
      cout << "Not familiar with the HDF4 type.. exiting.." << endl;
    }    
  }

  void * return_array(){
    switch( type_ ){
    case DFNT_INT8:
      return a_int8;
    case DFNT_UINT8:
      return a_uint8;
    case DFNT_INT16:
      return a_int16;
    case DFNT_UINT16:
      return a_uint16;
    case DFNT_INT32:
      return a_int32;
    case DFNT_UINT32:
      return a_uint32;
    case DFNT_FLOAT32:
      return a_float32;
    case DFNT_FLOAT64:
      return a_float64;
     default:
      cout << "Not familiar with the HDF4 type.. exiting.." << endl;
      return FAILURE;
    }  
  }

  double return_array_value(int i){

    //coverting to double. Not ideal.
    switch( type_ ){
    case DFNT_INT8:
      return a_int8[i];
    case DFNT_UINT8:
      return a_uint8[i];
    case DFNT_INT16:
      return a_int16[i];
    case DFNT_UINT16:
      return a_uint16[i];
    case DFNT_INT32:
      return a_int32[i];
    case DFNT_UINT32:
      return a_uint32[i];
    case DFNT_FLOAT32:
      return a_float32[i];
    case DFNT_FLOAT64:
      return a_float64[i];
    default:
      cout << "Not familiar with the HDF4 type.. exiting.." << endl;
      return FAILURE;
    }          
  }


  
};

/***************************************************************/
int read_hdf4(string file_name, Double_2D & data, char * data_name){

  //open the file
  int32 sd_id = SDstart(file_name.c_str(), DFACC_READ);
  if (sd_id == FAIL){
    cout << "Failed to open file:"<<file_name<< endl;
    return FAILURE;
  }
  
   //read the data block in the file
   int32 sds_index;
   sds_index = SDselect(sd_id,SDnametoindex(sd_id, data_name));
   if (sd_id == FAIL){
     cout << "Failed to find data block in file:"<<file_name<< endl;
     return FAILURE;
   }
   
   char sds_name[60];
   int32 rank, data_type, n_attrs;
   int32 dim_sizes[2];

   int status = SDgetinfo(sds_index, sds_name, &rank, 
			   dim_sizes, &data_type, &n_attrs);

   if(status){
     cout << "Failed getting data info"<< endl;
     return FAILURE;
   }
   if(rank!=2){
     cout << "Data block is not 2-dimensional"<< endl;
     return FAILURE;
   }

   hdf_anonymous_array array(data_type,dim_sizes[0]*dim_sizes[1] );

   int32 start[2] = {0};
   status = SDreaddata(sds_index, start, NULL, 
		       dim_sizes, array.return_array()  );
   if (status == FAIL){
     cout << "Could not get data block"<< endl;
     return FAILURE;
   }

   //close the SD readers
   status = SDendaccess(sds_index);
   status = SDend(sd_id);  


   //Fill return parameters
   if(data.get_size_x()==0)
     data.allocate_memory(dim_sizes[0],dim_sizes[1]);

   for(int i=0; i< dim_sizes[0]; ++i){
     for(int j=0; j< dim_sizes[1]; ++j){
       data.set(i,j,array.return_array_value(dim_sizes[1]*i+j));
     }
   }
   

   return SUCCESS;

 }
 
