#include <mpi.h>
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>

#include "particle.h"

/*The commonly refered leapfrog method is a type of predictor-corrector method.
we use the name predictor and corrector instead of kick and drift*/

/*predictor evovled to the i+1/2 time step from the i-1/2 step,
used to evolve velocity && (density || internal energy) (depending on the formulation)*/
void predictor(double NumoPar, double TimStp, std::vector<Particle *> VecoPtr2Par)
{
  std::cout << "mark 5.1 predictor start"<< std::endl;
  for (size_t i = 0; i < NumoPar; i++) {
    // ParIntEng += 0.5 * (TimStp[0]+TimStp[1]) * ParIntEngVel;
    // ParVel += 0.5 * (TimStp[0]+TimStp[1]) * ParAcc;

    VecoPtr2Par[i] -> AddIntEng(TimStp * (VecoPtr2Par[i] -> GetIntEngVel()));
    std::vector<double> temp_increment;
    temp_increment.clear();
    for (size_t j = 0; j < VecoPtr2Par[i] -> GetAcc().size(); j++) {
      temp_increment.push_back(TimStp * VecoPtr2Par[i] -> GetAcc()[j]);
    }
    VecoPtr2Par[i] -> AddVel(temp_increment);
    if(i % 10000 == 0){
    std::cout << "mark 5.1.1 this particle's timstp inteng acc vel: "<<VecoPtr2Par[i] ->GetIntEng()<<" ";
    for (auto acc : VecoPtr2Par[i] -> GetAcc()) {std::cout << acc << " ";}
    for (auto inc : temp_increment) {std::cout << inc << " ";}
    std::cout << std::endl;
    }
  }
}

/*corrector evovled to the i+1 time step from the i step,
used to evolve position && (density rate || internal energy rate) (depending on the formulation)*/
void corrector(double NumoNgb, double NumoPar, std::vector<Particle *> VecoPtr2Par, double &TimStp,
                 double Gam, double &SmtLenMax, double MshCof, std::vector<double> Bdr,
                 std::map<int, std::map<int, std::vector<Particle *>>> &Msh2Ptr2Par)
{
  std::vector<double> temp_increment;
  std::cout << "mark 4.1 corrector start: updating positions"<< std::endl;
  for (size_t i = 0; i < NumoPar; i++) {
    temp_increment.clear(); //initialise for each particle
    for (size_t j = 0; j < VecoPtr2Par[i] -> GetAcc().size(); j++) {
      temp_increment.push_back(TimStp * (VecoPtr2Par[i] -> GetAcc()[j]));
    }
    VecoPtr2Par[i] -> AddCor(temp_increment);
    VecoPtr2Par[i] -> ChkPrdBdr(Bdr);
    if(i % 10000 == 0){
    std::cout << "mark 4.1.1 this particle's cor inc: ";
    for (auto inc : temp_increment) {std::cout << inc << " ";}
    std::cout << std::endl;
    }
  }
  std::cout << "mark 4.2 generating meshes"<< std::endl;
  //calculate the smoothing length (search through neigbours again) and the Kernel weighted density
  generate_2dmesh(SmtLenMax, MshCof, NumoPar, VecoPtr2Par, Msh2Ptr2Par);
  std::cout << "mark 4.3 searching in meshes"<< std::endl;
  search_in_2dmesh(NumoNgb, NumoPar, VecoPtr2Par, Bdr, Msh2Ptr2Par, SmtLenMax, MshCof);
  std::cout << "mark 4.4 calculate smoothinglength using neighbours"<< std::endl;
  calculate_sml_den_prs(NumoNgb, NumoPar, VecoPtr2Par, Gam, SmtLenMax);
  std::cout << "mark 4.5 calculating acc intengvel using smoothinglength"<< std::endl;

  //calculate the acceleration at i+1 time step using coordinates pressure velocity etc.
  //calculate the internal energy rate at i+1 time step using coordinates pressure velocity etc.
  //calculate next time step
  calculate_acc_intengvel_timstp(NumoNgb, NumoPar, VecoPtr2Par, TimStp, Bdr);
  std::cout << "mark 4.6 corrector finished"<< std::endl;
}
