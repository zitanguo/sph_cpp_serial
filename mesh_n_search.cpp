#include <mpi.h>
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>

#include "particle.h"

// assumed the domain is from 0 to max
// return the map mapping from the lower left coners of the meshes to a vector of PIDs inside this mesh
void generate_2dmesh(double SmtLenMax, double MshCof, int NumoPar, std::vector<Particle *> VecoPtr2Par,
                     std::map<int, std::map<int, std::vector<Particle *>>> &Msh2Ptr2Par)
{
  Msh2Ptr2Par.clear(); //clear the map before using it
  int temp_x_idx, temp_y_idx; //temperary x and y indices
  //loop over particles and fill them into meshes
  for (int i = 0; i < NumoPar; i++)
  {
    temp_x_idx = floor(VecoPtr2Par[i] -> GetCor()[0] / (SmtLenMax * MshCof));
    temp_y_idx = floor(VecoPtr2Par[i] -> GetCor()[1] / (SmtLenMax * MshCof));
    Msh2Ptr2Par[temp_x_idx][temp_y_idx].push_back(VecoPtr2Par[i]); //add the pointer to the particle object to the dictionary
    if(i % 10000 == 0){
    std::cout << "mark 4.2.1 mesh x cor, coef, x idx: "<<VecoPtr2Par[i] -> GetCor()[0]
              <<" "<<(SmtLenMax * MshCof) <<" "<<temp_x_idx << std::endl;}
  }

}

void search_in_2dmesh(int NumoNgb, int NumoPar, std::vector<Particle *> VecoPtr2Par, std::vector<double> Bdr,
                      std::map<int, std::map<int, std::vector<Particle *>>> Msh2Ptr2Par, double SmtLenMax, double MshCof)
{
  double temp_dx = 0, temp_dy = 0, temp_distance = 0;
  double LenX = Bdr[2] - Bdr[0], LenY = Bdr[3] - Bdr[1];
  //forcely loop over every mesh
  size_t temp_i_total = ceil(LenX / (SmtLenMax * MshCof)),
         temp_j_total = ceil(LenY / (SmtLenMax * MshCof));
  std::cout << "mark 4.3.1 searching in "<<temp_i_total<<" * "<<temp_j_total<< " meshes "<<std::endl;
  //looping over every mesh
  for (size_t i = 0; i < temp_i_total; i++) {
    for (size_t j = 0; j < temp_j_total; j++) {
      //std::cout << "mark 2.2 size of this mesh "<<Msh2Ptr2Par[i][j].size() << std::endl;
      //looping over every particle in a mesh
      for (size_t n = 0; n < Msh2Ptr2Par[i][j].size(); n++) {
        //particle pairs in this mesh, don't need to check boundary condition
        for (size_t m = n + 1; m < Msh2Ptr2Par[i][j].size(); m++) {
          temp_distance = pow(Msh2Ptr2Par[i][j][n] -> GetCor()[0] - Msh2Ptr2Par[i][j][m] -> GetCor()[0], 2) +
                          pow(Msh2Ptr2Par[i][j][n] -> GetCor()[1] - Msh2Ptr2Par[i][j][m] -> GetCor()[1], 2);
          temp_distance = pow(temp_distance, 0.5);
          Msh2Ptr2Par[i][j][n] -> AddDisMap(temp_distance, Msh2Ptr2Par[i][j][m]);
          Msh2Ptr2Par[i][j][m] -> AddDisMap(temp_distance, Msh2Ptr2Par[i][j][n]);
        }


        //check the neigbouring meshes, periodic BC assumed
        size_t temp_i_mns_one = i - 1, temp_j_mns_one = j - 1,
               temp_i_pls_one = i + 1, temp_j_pls_one = j + 1;

        if (i == 0) {temp_i_mns_one = temp_i_total - 1;} // the left most bdr, i - 1 to right most
        else if (i == temp_i_total - 1) {temp_i_pls_one = 0;} // the right most bdr, i + 1 to the left most
        //assumed total j size is constant w. r. t. i
        if (j == 0) {temp_j_mns_one = temp_j_total - 1;} // the left most bdr, i - 1 to right most
        else if (j == temp_j_total - 1) {temp_j_pls_one = 0;} // the right most bdr, i + 1 to the left most


        //mth particle in the left mesh, i-1, j
        for (size_t m = 0; m < Msh2Ptr2Par[temp_i_mns_one][j].size(); m++)
        {
          temp_dx = std::abs(Msh2Ptr2Par[i][j][n] -> GetCor()[0] - Msh2Ptr2Par[temp_i_mns_one][j][m] -> GetCor()[0]);
          temp_dx = std::min(temp_dx, LenX - temp_dx);

          temp_dy = std::abs(Msh2Ptr2Par[i][j][n] -> GetCor()[1] - Msh2Ptr2Par[temp_i_mns_one][j][m] -> GetCor()[1]);
          temp_dy = std::min(temp_dy, LenY - temp_dy);

          temp_distance = sqrt(temp_dx * temp_dx + temp_dy * temp_dy);
          Msh2Ptr2Par[i][j][n] -> AddDisMap(temp_distance, Msh2Ptr2Par[temp_i_mns_one][j][m]);
          Msh2Ptr2Par[temp_i_mns_one][j][m] -> AddDisMap(temp_distance, Msh2Ptr2Par[i][j][n]);
        }

        //mth particle in the upper left mesh, i-1, j+1
        for (size_t m = 0; m < Msh2Ptr2Par[temp_i_mns_one][temp_j_pls_one].size(); m++)
        {
          temp_dx = std::abs(Msh2Ptr2Par[i][j][n] -> GetCor()[0] - Msh2Ptr2Par[temp_i_mns_one][temp_j_pls_one][m] -> GetCor()[0]);
          temp_dx = std::min(temp_dx, LenX - temp_dx);

          temp_dy = std::abs(Msh2Ptr2Par[i][j][n] -> GetCor()[1] - Msh2Ptr2Par[temp_i_mns_one][temp_j_pls_one][m] -> GetCor()[1]);
          temp_dy = std::min(temp_dy, LenY - temp_dy);

          temp_distance = sqrt(temp_dx * temp_dx + temp_dy * temp_dy);
          Msh2Ptr2Par[i][j][n] -> AddDisMap(temp_distance, Msh2Ptr2Par[temp_i_mns_one][temp_j_pls_one][m]);
          Msh2Ptr2Par[temp_i_mns_one][temp_j_pls_one][m] -> AddDisMap(temp_distance, Msh2Ptr2Par[i][j][n]);
        }

        //mth paritcle in the upper mesh, i, j+1
        for (size_t m = 0; m < Msh2Ptr2Par[i][temp_j_pls_one].size(); m++)
        {
          temp_dx = std::abs(Msh2Ptr2Par[i][j][n] -> GetCor()[0] - Msh2Ptr2Par[i][temp_j_pls_one][m] -> GetCor()[0]);
          temp_dx = std::min(temp_dx, LenX - temp_dx);

          temp_dy = std::abs(Msh2Ptr2Par[i][j][n] -> GetCor()[1] - Msh2Ptr2Par[i][temp_j_pls_one][m] -> GetCor()[1]);
          temp_dy = std::min(temp_dy, LenY - temp_dy);

          temp_distance = sqrt(temp_dx * temp_dx + temp_dy * temp_dy);
          Msh2Ptr2Par[i][j][n] -> AddDisMap(temp_distance, Msh2Ptr2Par[i][temp_j_pls_one][m]);
          Msh2Ptr2Par[i][temp_j_pls_one][m] -> AddDisMap(temp_distance, Msh2Ptr2Par[i][j][n]);
        }

        //mth paritcle in the upper mesh, i+1, j+1
        for (size_t m = 0; m < Msh2Ptr2Par[temp_i_pls_one][temp_j_pls_one].size(); m++)
        {
          temp_dx = std::abs(Msh2Ptr2Par[i][j][n] -> GetCor()[0] - Msh2Ptr2Par[temp_i_pls_one][temp_j_pls_one][m] -> GetCor()[0]);
          temp_dx = std::min(temp_dx, LenX - temp_dx);

          temp_dy = std::abs(Msh2Ptr2Par[i][j][n] -> GetCor()[1] - Msh2Ptr2Par[temp_i_pls_one][temp_j_pls_one][m] -> GetCor()[1]);
          temp_dy = std::min(temp_dy, LenY - temp_dy);

          temp_distance = pow(Msh2Ptr2Par[i][j][n] -> GetCor()[0] - Msh2Ptr2Par[temp_i_pls_one][temp_j_pls_one][m] -> GetCor()[0], 2) +
                          pow(Msh2Ptr2Par[i][j][n] -> GetCor()[1] - Msh2Ptr2Par[temp_i_pls_one][temp_j_pls_one][m] -> GetCor()[1], 2);
          temp_distance = sqrt(temp_dx * temp_dx + temp_dy * temp_dy);
          Msh2Ptr2Par[i][j][n] -> AddDisMap(temp_distance, Msh2Ptr2Par[temp_i_pls_one][temp_j_pls_one][m]);
          Msh2Ptr2Par[temp_i_pls_one][temp_j_pls_one][m] -> AddDisMap(temp_distance, Msh2Ptr2Par[i][j][n]);
        }
      }
    }
  }
}
