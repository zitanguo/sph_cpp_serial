#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>

#include "particle.h"

double kernel_2d3(double r, double h){
  if(r >= (double) 2 * h ) {return 0.0;}
  else{
    double r_over_h = r / h;
    double norm = 10.0/(7.0 * M_PI * h * h);
    return ((r_over_h > 1) ? norm / 4.0 * pow((2 - r_over_h), 3) :
            norm * (1 - 1.5 * pow(r_over_h, 2) + 0.75 * pow(r_over_h, 3)));
  }
}

double kernel_2d3_derivative_r(double r, double h){
  if(r >= (double) 2 * h) {return 0.0;}
  else{
    double r_over_h = r / h;
    double norm = 10.0/(7.0 * M_PI * h * h * h);
    return ((r_over_h > 1) ? -3 * norm / 4.0 * pow((2 - r_over_h), 2) :
            norm * (- 3.0 * r_over_h + 9.0/4.0 * pow(r_over_h, 2)));
  }
}

//use gnuplot to have a quick check on the kernel shapes
/*
int main() {

  std::ofstream myfile;
  myfile.open ("kernel.txt");
  double dr = 1./1000., h = 0.5;
  for (size_t i = 0; i < 1000; i++) {
    myfile <<((double)i * dr)<<" "<<kernel_2d3((double)i * dr, h) <<" "<<kernel_2d3_derivative_r((double)i * dr, h)<<"\n";
  }
  myfile.close();
  return 0;
}
*/
