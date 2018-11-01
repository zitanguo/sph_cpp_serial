#include <mpi.h>
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <iterator>

#include "particle.h"

void calculate_sml_den_prs(int NumoNgb, int NumoPar, std::vector<Particle *> VecoPtr2Par,
                           double gamma, double &SmtLenMax)
{
  for (size_t i = 0; i < NumoPar; i++)
  {
    //*****find the smoothing length first*****
    //cpp maps are sorted, so just use the # of ngbth's distance as the smt
    std::map<double, Particle *>::iterator temp_stopper = VecoPtr2Par[i] -> GetDisMap().begin();
    if(NumoNgb <= VecoPtr2Par[i] -> GetDisMap().size())
    {
      std::advance(temp_stopper, NumoNgb - 1); //move the iterator to the ngb th distance}
      VecoPtr2Par[i] -> SetSmtLen(temp_stopper -> first);  //"fisrt" is the key or the distance, set to the smt
    }
    else
    {
      std::cout << "warning: insufficient neighbours: " << VecoPtr2Par[i] -> GetDisMap().size()<<std::endl;
      temp_stopper = VecoPtr2Par[i] -> GetDisMap().end();
      VecoPtr2Par[i] -> SetSmtLen(temp_stopper -> first);  //"fisrt" is the key or the distance, set to the smt
    }
    if(temp_stopper -> first > SmtLenMax){SmtLenMax = temp_stopper -> first;}


    //*****then find the density using the just found smoothing length*****
    int temp_counter = 0; //count the times of iteration
    double temp_density = 0.0;
    //this for loop needs -std=c++17
    for(auto const& [Dis, Ptr2Par] : VecoPtr2Par[i] -> GetDisMap()){
      temp_counter++;
      temp_density += Ptr2Par -> GetMas() * kernel_2d3(Dis, temp_stopper -> first);
      if (temp_counter >= NumoNgb) break;
    }
    VecoPtr2Par[i] -> SetDen(temp_density);

    //*****then set the pressure using density and internal energy*****
    VecoPtr2Par[i] -> CalPrs(gamma);
    if(i % 10000 == 0){
    std::cout << "mark 4.4.1 this particle's h, rho, prs: "<<temp_stopper -> first <<" "<<temp_density
              <<" "<<VecoPtr2Par[i] -> GetPrs() <<std::endl;}
  }
}

//This loop is dependent on the pressure, density and sml of every particle so has to be carried out after the first one.
void calculate_acc_intengvel_timstp(int NumoNgb, int NumoPar, std::vector<Particle *> VecoPtr2Par, double &TimStp, std::vector<double> Bdr)
{

  for (size_t i = 0; i < NumoPar; i++)
  {
    VecoPtr2Par[i] -> CalAccIntEngVel(NumoNgb, Bdr);
    //clear the map here because the assignment isn't particlely linear (it's meshly linear)
    VecoPtr2Par[i] -> GetDisMap().clear();
    if(i % 10000 == 0){
    std::cout << "mark 4.5.1 this particle's engvel, acc: "<<VecoPtr2Par[i] -> GetIntEngVel()<<" ";
    for (auto acc : VecoPtr2Par[i] -> GetAcc()) {std::cout << acc << " ";}
    std::cout << std::endl;
    }

    //calculate timestep
    double temp_time_step = pow(VecoPtr2Par[i] -> GetSmtLen() / VecoPtr2Par[i] -> GetAccAbs(), 0.5);
    // std::cout <<i <<" "<<VecoPtr2Par[i] -> GetDen()<<" "<<VecoPtr2Par[i] -> GetMas() <<" "<<VecoPtr2Par[i] -> GetIntEng()<<" "
    //           <<VecoPtr2Par[i] -> GetSmtLen()<<" "<<VecoPtr2Par[i] -> GetPrs()<<" "<<VecoPtr2Par[i] -> GetAcc()[0]<<" "
    //           <<VecoPtr2Par[i] -> GetAcc()[1]<<" "<< VecoPtr2Par[i] -> GetAccAbs()<<std::endl;
    //std::cout <<temp_time_step<<" " <<TimStp<<std::endl;
    if( temp_time_step < TimStp)
      TimStp = temp_time_step;
  }
  std::cout << "mark 4.5.2 finished calculating acc with timestep "<<TimStp<<std::endl;
}
