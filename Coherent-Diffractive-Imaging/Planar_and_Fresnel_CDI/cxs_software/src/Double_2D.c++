#include "string.h"
#include "Double_2D.h"
  
Double_2D::Double_2D():nx(0),ny(0){
}

Double_2D::Double_2D(int x_size, int y_size){
    allocate_memory(x_size,y_size);
}

// Destructor
Double_2D::~Double_2D(){
  if(nx > 0 ){
    delete [] array;
  }  
}

//I don't think making a 2-D array 
//from a 1-D array actually saves much CPU
//but at least copying might be faster.
void Double_2D::allocate_memory(int x_size, int y_size){
  nx = x_size;
  ny = y_size;
  array = new double[nx*ny];
  for(int i=0; i<nx; i++)
    for(int j=0; j<ny; j++)
      array[i*ny+j]=0;
}


//fast copy, no error checking
void Double_2D::copy(const Double_2D & double_array){
  memcpy(array,double_array.array, sizeof(double)*nx*ny);
}


double Double_2D::get_sum() const {
  double total = 0;
  for(int i=0; i<nx; i++)
    for(int j=0; j<ny; j++)
      total+=array[i*ny+j];
  return total;
};

double Double_2D::get_max() const {
  if(nx==0||ny==0)
    return 0;

  double max = array[0];
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      if(array[i*ny+j]>max)
	max = array[i*ny+j];
    }
  }
  return max;
};


double Double_2D::get_min() const {
  if(nx==0||ny==0)
    return 0;

  double min = array[0];
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      if(array[i*ny+j]<min)
	min = array[i*ny+j];
    }
  }
  return min;
};
