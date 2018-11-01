#include <mpi.h>
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>

#include "particle.h"

/*----------------------------------------------------------------------------
 * A hard coded programme to read initial H5 file. H5 doesn't support cpp so use c-style arrays here
 */


int read_initial_h5(char* file_name, char* group_name,
                     float* &ParCor, float* &ParDen, float* &ParIntEng,
                     float* &ParMas, int* &ParIds, float* &ParSmtLen, float* &ParVel)
{
  size_t rank, ndims, num_par, num_elem;
  const char * dataset_name_list[] = {"Coordinates", "Density", "InternalEnergy", "Masses", "ParticleIDs", "SmoothingLength", "Velocities"};
  //printf("%d\n", sizeof(dataset_name_list)/sizeof(*dataset_name_list));
  hid_t   file_id, group_id, dataset_id, dataspace_id, memspace_id, datatype_id;
  hsize_t* dims;



  //Open the specified h5 file, group, dataset and dataspace
  file_id = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
  group_id = H5Gopen2(file_id, group_name, H5P_DEFAULT);
  for (size_t i = 0; i < sizeof(dataset_name_list)/sizeof(*dataset_name_list); i++)
  {
    printf("\n reading %s\n", dataset_name_list[i]);
    dataset_id = H5Dopen2(group_id, dataset_name_list[i], H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    datatype_id = H5Dget_type(dataset_id);


    rank = H5Sget_simple_extent_ndims(dataspace_id);
    //dims = (hsize_t*) malloc(rank * sizeof(hsize_t));
    dims = new hsize_t[rank];
    ndims = H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
    memspace_id = H5Screate_simple(rank, dims, NULL);

    num_elem = H5Sget_simple_extent_npoints(dataspace_id);

    for (size_t i = 0; i < rank; i++) {
      printf("dim %ld: %ld ", i, dims[i]);
    }
    printf("total: %ld\n", num_elem);

    //A hard coded if conditions to identify different type of variables
    if(dataset_name_list[i] == "Coordinates"){
      ParCor = new float[num_elem];
      H5Dread(dataset_id, datatype_id, memspace_id, dataspace_id, H5P_DEFAULT, ParCor);
    }
    else if(dataset_name_list[i] == "Density"){
      //double ParDen[dims[0]];
      ParDen = new float[num_elem];
      H5Dread(dataset_id, datatype_id, memspace_id, dataspace_id, H5P_DEFAULT, ParDen);
    }
    else if(dataset_name_list[i] == "InternalEnergy"){
      //double ParIntEng[dims[0]];
      ParIntEng = new float[num_elem];
      H5Dread(dataset_id, datatype_id, memspace_id, dataspace_id, H5P_DEFAULT, ParIntEng);
    }
    else if(dataset_name_list[i] == "Masses"){
      //double ParMas[dims[0]];
      ParMas = new float[num_elem];
      H5Dread(dataset_id, datatype_id, memspace_id, dataspace_id, H5P_DEFAULT, ParMas);
    }
    else if(dataset_name_list[i] == "ParticleIDs"){
      //double ParIds[dims[0]];
      ParIds = new int[num_elem];
      H5Dread(dataset_id, datatype_id, memspace_id, dataspace_id, H5P_DEFAULT, ParIds);
    }
    else if(dataset_name_list[i] == "SmoothingLength"){
      //double ParSmtLen[dims[0]];
      ParSmtLen = new float[num_elem];
      H5Dread(dataset_id, datatype_id, memspace_id, dataspace_id, H5P_DEFAULT, ParSmtLen);
    }
    else if(dataset_name_list[i] == "Velocities"){
      //double ParVel[dims[0]][dims[1]];
      ParVel = new float[num_elem];
      H5Dread(dataset_id, datatype_id, memspace_id, dataspace_id, H5P_DEFAULT, ParVel);
    }
    else
      printf("Unexpected rank %d in group %s dataset %s\n", rank, group_name, dataset_name_list[i]);

    //ptr_out = (double *)malloc((int)num_elem * sizeof(double));
    //H5Dread(dataset_id, datatype_id, memspace_id, dataspace_id, H5P_DEFAULT, ptr_out);

    H5Tclose(datatype_id);
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    num_par = dims[0];
    delete[] dims;
  }
  H5Dclose(dataset_id);
  H5Gclose(group_id);
  H5Fclose(file_id);
  std::cout << "Num of Par is: " <<num_par<<'\n';
  return num_par;
}

//***** write to h5 files*****
void write_h5(char* file_name, std::vector<Particle *> VecoPtr2Par)

{
  int rank, ndims, num_par;
  const char * dataset_name_list[] = {"Coordinates", "Density", "InternalEnergy", "Masses", "ParticleIDs", "SmoothingLength", "Velocities"};
  //printf("%d\n", sizeof(dataset_name_list)/sizeof(*dataset_name_list));
  hid_t   file_id, group_id, dataset_id, dataspace_id, memspace_id, datatype_id;
  hsize_t num_elem;
  hsize_t* dims;

  int dim0 = VecoPtr2Par.size(), dim1 = VecoPtr2Par[0] -> GetCor().size();
  double ParDen[VecoPtr2Par.size()], ParIntEng[VecoPtr2Par.size()], ParMas[VecoPtr2Par.size()],
         ParIds[VecoPtr2Par.size()], ParSmtLen[VecoPtr2Par.size()];
  double ParCor[VecoPtr2Par.size()][VecoPtr2Par[0] -> GetCor().size()],
         ParVel[VecoPtr2Par.size()][VecoPtr2Par[0] -> GetVel().size()];

  for (size_t i = 0; i < VecoPtr2Par.size(); i++) {
    ParDen[i] = VecoPtr2Par[i] -> GetDen();
    ParIntEng[i] = VecoPtr2Par[i] -> GetIntEng();
    ParMas[i] = VecoPtr2Par[i] -> GetMas();
    ParIds[i] = VecoPtr2Par[i] -> GetIds();
    ParSmtLen[i] = VecoPtr2Par[i] -> GetSmtLen();
    for (size_t j = 0; j < VecoPtr2Par[0] -> GetCor().size(); j++) {
      ParCor[i][j] = VecoPtr2Par[i] -> GetCor()[j];
      ParVel[i][j] = VecoPtr2Par[i] -> GetVel()[j];
    }
  }



  //Open the specified h5 file, group, dataset and dataspace
  file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  //group_id = H5Gopen2(file_id, group_name, H5P_DEFAULT);
  for (size_t i = 0; i < sizeof(dataset_name_list)/sizeof(*dataset_name_list); i++)
  {
    //A hard coded if conditions to identify different type of variables
    if(dataset_name_list[i] == "Coordinates"){
      rank = 2;
      dims = new hsize_t[rank];
      dims[0] = dim0;
      dims[1] = dim1;
      dataspace_id = H5Screate_simple(rank, dims, NULL);
      datatype_id = H5Tcopy(H5T_NATIVE_DOUBLE);
      H5Tset_order(datatype_id, H5T_ORDER_LE);
      dataset_id = H5Dcreate2(file_id, dataset_name_list[i], datatype_id, dataspace_id,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, ParCor);
      H5Sclose(dataspace_id);
      H5Tclose(datatype_id);
      delete[] dims;
      H5Dclose(dataset_id);

    }
    else if(dataset_name_list[i] == "Density"){
      rank = 1;
      dims = new hsize_t[rank];
      dims[0] = dim0;
      //dims[1] = dim1;
      dataspace_id = H5Screate_simple(rank, dims, NULL);
      datatype_id = H5Tcopy(H5T_NATIVE_DOUBLE);
      H5Tset_order(datatype_id, H5T_ORDER_LE);
      dataset_id = H5Dcreate2(file_id, dataset_name_list[i], datatype_id, dataspace_id,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, ParDen);
      H5Sclose(dataspace_id);
      H5Tclose(datatype_id);
      delete[] dims;
      H5Dclose(dataset_id);
    }
    else if(dataset_name_list[i] == "InternalEnergy"){
      rank = 1;
      dims = new hsize_t[rank];
      dims[0] = dim0;
      //dims[1] = dim1;
      dataspace_id = H5Screate_simple(rank, dims, NULL);
      datatype_id = H5Tcopy(H5T_NATIVE_DOUBLE);
      H5Tset_order(datatype_id, H5T_ORDER_LE);
      dataset_id = H5Dcreate2(file_id, dataset_name_list[i], datatype_id, dataspace_id,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, ParIntEng);
      H5Sclose(dataspace_id);
      H5Tclose(datatype_id);
      delete[] dims;
      H5Dclose(dataset_id);
    }
    else if(dataset_name_list[i] == "Masses"){
      rank = 1;
      dims = new hsize_t[rank];
      dims[0] = dim0;
      //dims[1] = dim1;
      dataspace_id = H5Screate_simple(rank, dims, NULL);
      datatype_id = H5Tcopy(H5T_NATIVE_DOUBLE);
      H5Tset_order(datatype_id, H5T_ORDER_LE);
      dataset_id = H5Dcreate2(file_id, dataset_name_list[i], datatype_id, dataspace_id,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, ParMas);
      H5Sclose(dataspace_id);
      H5Tclose(datatype_id);
      delete[] dims;
      H5Dclose(dataset_id);
    }
    else if(dataset_name_list[i] == "ParticleIDs"){
      rank = 1;
      dims = new hsize_t[rank];
      dims[0] = dim0;
      //dims[1] = dim1;
      dataspace_id = H5Screate_simple(rank, dims, NULL);
      datatype_id = H5Tcopy(H5T_NATIVE_INT);
      H5Tset_order(datatype_id, H5T_ORDER_LE);
      dataset_id = H5Dcreate2(file_id, dataset_name_list[i], datatype_id, dataspace_id,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, ParIds);
      H5Sclose(dataspace_id);
      H5Tclose(datatype_id);
      delete[] dims;
      H5Dclose(dataset_id);
    }
    else if(dataset_name_list[i] == "SmoothingLength"){
      rank = 1;
      dims = new hsize_t[rank];
      dims[0] = dim0;
      //dims[1] = dim1;
      dataspace_id = H5Screate_simple(rank, dims, NULL);
      datatype_id = H5Tcopy(H5T_NATIVE_DOUBLE);
      H5Tset_order(datatype_id, H5T_ORDER_LE);
      dataset_id = H5Dcreate2(file_id, dataset_name_list[i], datatype_id, dataspace_id,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, ParSmtLen);
      H5Sclose(dataspace_id);
      H5Tclose(datatype_id);
      delete[] dims;
      H5Dclose(dataset_id);
    }
    else if(dataset_name_list[i] == "Velocities"){
      rank = 2;
      dims = new hsize_t[rank];
      dims[0] = dim0;
      dims[1] = dim1;
      dataspace_id = H5Screate_simple(rank, dims, NULL);
      datatype_id = H5Tcopy(H5T_NATIVE_DOUBLE);
      H5Tset_order(datatype_id, H5T_ORDER_LE);
      dataset_id = H5Dcreate2(file_id, dataset_name_list[i], datatype_id, dataspace_id,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, ParVel);
      H5Sclose(dataspace_id);
      H5Tclose(datatype_id);
      delete[] dims;
      H5Dclose(dataset_id);
    }


  }
  H5Fclose(file_id);
  //std::cout << "Num of Par is: " <<num_par<<'\n';
}
