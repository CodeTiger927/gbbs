#include "ligra/pbbslib/dyn_arr.h"


//No need for vertex_data because it would store dyn_arr of neighbors anyways
struct dyn_arr_sym_vertex {
  dyn_arr<uintE> neighbors; 
  uintE deg;
  //Doesn't care about edges because it has neighbors and weight doesn't matter

  void add_neighbor(uintE v) {
    deg++;
    
  }

  void remove_neighbor(uintE v) {
    //Change entry for v into last neighbor of vertex
    //Set last entry to empty value
  }
};


//Graph is unweighted
struct dyn_arr_sym_graph {
  using edge_type = uintE; //No need for tuple with weight

  dyn_arr<dyn_arr_sym_vertex> vertices;

  size_t n; //Number of vertices
  size_t m; //Number of edges

  //Called to delete graph
  std::function<void()> _deletion_fn;


  dyn_arr_sym_graph(std::function<void()> _deletion_fn):
    deletion_fn(_deletion_fn),
    n(0), //Graph starts with no vertices or edges
    m(0),
    vertices(dyn_arr(20)) //What is the best size
  {}



  /* Why is this needed but only for edges
  #ifndef TWOSOCKETNVM
    dyn_arr_sym_vertex get_vertex(uintE i) {
    
    }
  #else
    dyn_arr_sym_vertex get_vertex(uintE i) {

    }
  #endif
  */




  dyn_arr_sym_vertex get_vertex(uintE i) {
    return vertices.A[i];
    
  }

  void del() {
    deletion_fn();
  }

  void add_vertex(uintE i) {

    vertices.push_back(i);
  }


  void remove_vertex(uintE i) {

    //Maybe just set i to empty value because the indices are needed to maintain ID
  }

  void add_edge(uintE u, uintE v) {


  }


  void remove_edge(uintE u, uintE v) {
    
    
    
  }



    


  





};
