#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <vector>
#include <numeric>
#include <map>
#include <iostream>
#include <iterator>
#include <math.h>
#include <algorithm>

#include "particle.h"


//Particle constructor
//Particle::Particle(double Rx, double Ry, double Rz, double Density, double InternalEnergy, double Masses,
//                   int ParticleIDs, double SmoothingLength, double Vx, double Vy, double Vz)

Particle::Particle(double Rx, double Ry, double Density, double InternalEnergy, double Masses,
                   int ParticleIDs, double SmoothingLength, double Vx, double Vy, double gamma)
{
  //Cor = {Rx, Ry, Rz};
  Cor = {Rx, Ry};
  Mas = Masses;
  IntEng = InternalEnergy;
  SmtLen = SmoothingLength;
  Ids = ParticleIDs;
  Den = Density;
  //Vel = {Vx, Vy, Vz};
  Vel = {Vx, Vy};
  Prs = (gamma - 1) * Den * IntEng;
  Acc = {0.0, 0.0};
}

//retrieve private values
int Particle::GetIds() {return Ids;}
double Particle::GetDen() {return Den;}
double Particle::GetMas() {return Mas;}
double Particle::GetIntEng() {return IntEng;}
double Particle::GetSmtLen() {return SmtLen;}
double Particle::GetPrs() {return Prs;}
double Particle::GetIntEngVel() {return IntEngVel;}

double Particle::GetAccAbs()
{return pow(std::inner_product(Acc.begin(), Acc.end(), Acc.begin(), 0.0), 0.5);} //abs value of this vector

std::vector<double> Particle::GetCor(){return Cor;}
std::vector<double> Particle::GetVel(){return Vel;}
std::vector<double> Particle::GetAcc(){return Acc;}
std::map<double, Particle *> Particle::GetDisMap(){return Dis2Ptr2Par;}

//set private values
void Particle::SetDen(double den) {Den = den;}
void Particle::SetSmtLen(double len) {SmtLen = len;}

//modify values
void Particle::AddDisMap(double dis, Particle *par) {Dis2Ptr2Par[dis] = par;}
void Particle::AddIntEng(double increment) {IntEng += increment;}
void Particle::AddVel(std::vector<double> increment)
{
  for (size_t i = 0; i < Vel.size(); i++)
  {
    Vel[i] += increment[i];
  }
}
void Particle::AddCor(std::vector<double> increment)
{
  for (size_t i = 0; i < Cor.size(); i++)
  {
    Cor[i] += increment[i];
  }
}

void Particle::ChkPrdBdr(std::vector<double> Bdr) // Bdr arranged in Xmin, Ymin, Zmin (if avail.),Xmax, Ymax, Zmax (if avail.)
{
  for (size_t i = 0; i < Cor.size(); i++)
  {
    if(Cor[i] < Bdr[i]) Cor[i] += (Bdr[Cor.size()+i] - Bdr[i]); //left bdr
    else if (Cor[i] > Bdr[Cor.size()+i]) Cor[i] -= (Bdr[Cor.size()+i] - Bdr[i]); //right bdr
  }
}

//calculate derived values
void Particle::CalPrs(double gamma) {Prs = (gamma - 1) * Den * IntEng;} //adiabatic EoS
//calculate the acceleration and internal energy change rate in one loop
void Particle::CalAccIntEngVel(double NumoNgb, std::vector<double> Bdr)
{
  //clear the original values
  std::fill(Acc.begin(), Acc.end(), 0.0);
  IntEngVel = 0.0;

  int temp_counter = 0;
  //this loop needs -std=c++17
  for(auto const& [Dis, Ptr2Par] : Dis2Ptr2Par)
  {
    double temp_acc_magnitude = 0.0, temp_intengvel_magnitude = 0.0,
           temp_intengvel_product = 0.0, temp_di = 0.0; //temp_di is temp dx/dy/etc.
    temp_counter++;
    temp_acc_magnitude = - Ptr2Par -> GetMas() / Dis *
                         (Prs/Den/Den * kernel_2d3_derivative_r(Dis,SmtLen) +
                         Ptr2Par -> GetPrs()/Ptr2Par -> GetDen()/Ptr2Par -> GetDen() *
                         kernel_2d3_derivative_r(Dis, Ptr2Par -> GetSmtLen()));

    temp_intengvel_magnitude =  Ptr2Par -> GetMas() / Dis * kernel_2d3_derivative_r(Dis, SmtLen);

    //periodic BC assumed
    for (size_t i = 0; i < Acc.size(); i++)
    {
      //ri -rj with direction preserved, |ri -rj| < L - |ri-rj| ?  ri -rj : (ri-rj<0? ri-rj+L:ri-rj-L)
      temp_di = Cor[i] - Ptr2Par -> GetCor()[i];
      temp_di = std::abs(temp_di) < std::abs(Bdr[Acc.size()+i] - Bdr[i] - std::abs(temp_di)) ?
                temp_di : (temp_di < 0 ? temp_di + Bdr[Acc.size()+i] - Bdr[i]: temp_di - Bdr[Acc.size()+i] + Bdr[i]);

      Acc[i] += temp_di * temp_acc_magnitude;
      //vector product of (ri - rj) * (vi - vj)
      temp_intengvel_product += temp_di * (Vel[i] - Ptr2Par -> GetVel()[i]);
    }
    IntEngVel += temp_intengvel_magnitude * temp_intengvel_product;
    if (temp_counter >= NumoNgb) break;
  }
  IntEngVel *= Prs / Den / Den; // this coefficient is for the ith particle itself, indepent on it's jth neigbour
  //std::cout << "intengvel acc "<<IntEngVel<<" "<<Acc[0]<<" "<<Acc[1] <<std::endl;
}

//destructor
Particle::~Particle () {std::cout << "Destructor called" << std::endl;}
