#include <mpi.h>
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>

class   quad_node{
private:
  int level; // the level of this node, initialised to zero
  std::vector<Particle *> LocPtr2Par; //pointers to particles contained in this node
  std::vector<double> LocBdr; //Boundary values, in the order of xmin, ymin, xmax, ymax
  quad_node *NodNW;     // child node representing northwest quadrant
  quad_node *NodNE;     // child node representing northeast quadrant
  quad_node *NodSW;     // child node representing southwest quadrant
  quad_node *NodSE;     // child node representing southeast quadrant
  // class quad_node *Prt;        //parent node of the current node

public:
  quad_node(int level, std::vector<Particle *> LocPtr2Par, std::vector<double> LocBdr);
  bool quan_node_check_particles() {return (bool) LocParIds.size()};
  void quad_node_split(); //split to furhter nodes
  virtual  ~quad_node();
};

  quad_node()
  {
    xmin = left;
    xmax = right;
    ymin = down;
    ymax = up;
    LocParIds = id_vector;
  }

  void quad_node::quad_node_split(){

  }
