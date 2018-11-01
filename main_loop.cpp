#include <iostream>
#include <vector>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <ctime>

#include "particle.h"

int main(int argc, char const *argv[]) {
  clock_t begin = clock();

  int Dim = 2; // Default dim is two
  int IcfDim = 3; // The initial condition file dimension
  int DatIdx = 1; //data index
  size_t NumoPar;
  double Gam = 5./3.; //adiabatic gamma
  double TotTim = 0;
  double TimLim = 10;
  double TimStp = 0.002; //use a large number
  double DelTim = 0.1; //data dump time step
  double NumoNgb = 20;
  double Xmin = 0, Xmax = 1;
  double Ymin = 0, Ymax = 1;
  double SmtLenMax = 0, MshCof = 2.0;
  std::vector<double> Bdr = {Xmin, Ymin, Xmax, Ymax};


  //some how the h5 api supports the c-style arrays only, read cor vel into 1d arrays and change to 2d later
  //double *ParCor, *ParDen, *ParIntEng, *ParMas, *ParSmtLen, *ParVel;
  float *ParCor, *ParDen, *ParIntEng, *ParMas, *ParSmtLen, *ParVel;
  int* ParIds;
  char* file_name = "../kh_mcnally_2d_ics.hdf5";
  char* group_name = "PartType0";

  std::cout << "mark 1 reading the initial condition file" << std::endl;
  //reading the initial h5 files into c-style arrays
  NumoPar = read_initial_h5(file_name, group_name, ParCor, ParDen, ParIntEng, ParMas, ParIds, ParSmtLen, ParVel);
  std::vector<Particle *> VecoPtr2Par; //vector of pointers to particles

  std::cout << "mark 2 assigning initial values" << std::endl;
  for (size_t i = 0; i < NumoPar; i++) {
    VecoPtr2Par.push_back(new Particle(ParCor[IcfDim*i], ParCor[IcfDim*i+1],
                                    ParDen[i], ParIntEng[i], ParMas[i], ParIds[i], ParSmtLen[i],
                                    ParVel[IcfDim*i], ParVel[IcfDim*i+1], Gam));
    if(ParSmtLen[i] > SmtLenMax){SmtLenMax = ParSmtLen[i];}
  }
  delete[] ParCor, ParDen, ParIntEng, ParMas, ParIds, ParSmtLen, ParVel;
  // mesh to a vector of pointers to the particles inside this mesh
  std::map<int, std::map<int, std::vector<Particle *>>> Msh2Ptr2Par;

  std::cout << "mark 3 initialising variables with max h "<<SmtLenMax<< std::endl;
  //*****initialising for time evolution*****
  //need to construct the mesh and find neighbors first in order to calculate acceleration
  std::cout << "mark 3.1 generating meshes "<< std::endl;
  generate_2dmesh(SmtLenMax, MshCof, NumoPar, VecoPtr2Par, Msh2Ptr2Par);
  std::cout << "mark 3.2 searching in meshes "<< std::endl;
  search_in_2dmesh(NumoNgb, NumoPar, VecoPtr2Par, Bdr, Msh2Ptr2Par, SmtLenMax, MshCof);
  std::cout << "mark 3.3 calculating the initial acceleration with timestep "<<TimStp << std::endl;
  // then calculate the initial acc and internal energy change rate at this position
  calculate_acc_intengvel_timstp(NumoNgb, NumoPar, VecoPtr2Par, TimStp, Bdr);
  std::cout << "mark 3.4 initilly evolving velocity "<< std::endl;
  predictor(NumoPar, 0.5 * TimStp, VecoPtr2Par);
  std::cout << "mark 3.5 finished initialising with timestep "<<TimStp<< std::endl;
  TotTim += 0.5 * TimStp;

  do {

    corrector(NumoNgb, NumoPar, VecoPtr2Par, TimStp, Gam, SmtLenMax, MshCof, Bdr, Msh2Ptr2Par);

    predictor(NumoPar, TimStp, VecoPtr2Par);

    //add data dump here
    if(TotTim >= ((double) DatIdx * DelTim))
    {
      std::string out_file_name = "out_data_"+std::to_string(DatIdx)+".hdf5";
      std::cout << "data dumping to "<<out_file_name<< std::endl;
      char* out_file_cname = strdup(out_file_name.c_str());
      write_h5(out_file_cname, VecoPtr2Par);
      free(out_file_cname);
      DatIdx++;
    }
    clock_t now = clock();
    std::cout << "h max   simulation time   system time" << std::endl;
    std::cout <<SmtLenMax<<"   "<<TotTim <<"   "<<double(now - begin) / CLOCKS_PER_SEC<<std::endl;
    TotTim += TimStp;
  } while(TotTim < TimLim);

  clock_t end = clock();
  std::cout << "Finised using " <<double(end - begin) / CLOCKS_PER_SEC<<" s"<<std::endl;

  return 0;

}
