#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <vector>
#include <numeric>
#include <map>
#include <iostream>
#include <iterator>


class Particle {
private:

  double Den, Mas,  IntEng, SmtLen, Ids;
  std::vector<double> Cor, Vel, Acc;
  double Prs;
  double IntEngVel;
  std::map<double, Particle *> Dis2Ptr2Par;

public:
  //Constructor
  // Particle (double Rx, double Ry, double Rz, double Density, double InternalEnergy, double Masses,
  //           int ParticleIDs, double SmoothingLength, double Vx, double Vy, double Vz);

  Particle (double Rx, double Ry, double Density, double InternalEnergy, double Masses,
            int ParticleIDs, double SmoothingLength, double Vx, double Vy, double gamma);
  //retrieve private values
  int GetIds();
  double GetDen();
  double GetMas();
  double GetIntEng();
  double GetSmtLen();
  double GetPrs();
  double GetIntEngVel();
  double GetAccAbs();
  std::vector<double> GetCor();
  std::vector<double> GetVel();
  std::vector<double> GetAcc();
  std::map<double, Particle *> GetDisMap();
  //set private values
  void SetDen(double den);
  void SetSmtLen(double len);
  //modify values
  void AddDisMap(double dis, Particle *par);
  void AddIntEng(double increment);
  void AddVel(std::vector<double> increment);
  void AddCor(std::vector<double> increment);
  void ChkPrdBdr(std::vector<double> Bdr);
  //calculate derived values
  void CalPrs(double gamma);
  void CalAccIntEngVel(double NumoNgb, std::vector<double> Bdr);


  virtual ~Particle ();
};

int read_initial_h5(char* file_name, char* group_name,
                     float* &ParCor, float* &ParDen, float* &ParIntEng,
                     float* &ParMas, int* &ParIds, float* &ParSmtLen, float* &ParVel);
void write_h5(char* file_name, std::vector<Particle *> VecoPtr2Par);

double kernel_2d3(double r, double h);
double kernel_2d3_derivative_r(double r, double h);

void generate_2dmesh(double SmtLenMax, double MshCof, int NumoPar, std::vector<Particle *> VecoPtr2Par,
                     std::map<int, std::map<int, std::vector<Particle *>>> &Msh2Ptr2Par);

void search_in_2dmesh(int NumoNgb, int NumoPar, std::vector<Particle *> VecoPtr2Par, std::vector<double> Bdr,
                     std::map<int, std::map<int, std::vector<Particle *>>> Msh2Ptr2Par, double SmtLenMax, double MshCof);

void calculate_sml_den_prs(int NumoNgb, int NumoPar, std::vector<Particle *> VecoPtr2Par, double gamma, double &SmtLenMax);
void calculate_acc_intengvel_timstp(int NumoNgb, int NumoPar, std::vector<Particle *> VecoPtr2Par,
                                    double &TimStp, std::vector<double> Bdr);

void predictor(double NumoPar, double TimStp, std::vector<Particle *> VecoPtr2Par);
void corrector(double NumoNgb, double NumoPar, std::vector<Particle *> VecoPtr2Par, double &TimStp,
               double Gam, double &SmtLenMax, double MshCof, std::vector<double> Bdr,
               std::map<int, std::map<int, std::vector<Particle *>>> &Msh2Ptr2Par);
