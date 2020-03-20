#include <iostream>
#include <map>
#include <set>
#include <limits>
#include <vector>
#include <math.h>
#include "ligra/ligra.h"
#include "dynamic_symmetric_graph.h"
#include "dyn_arr_sym_graph.h"
#include "ligra/pbbslib/sparse_table.h"
#include "ligra/pbbslib/dyn_arr.h"

// New Code for HSet
// For hashing sparse table

//Already defined
//struct hash_uintE {
  //inline size_t operator () (const uintE & a) {return pbbs::hash64_2(a);}
//};

//TODO: same notes as on h-index2

struct HSet {

  symmetric_graph<symmetric_vertex, pbbs::empty>* G; //Graph

  sparse_table<uintE, pbbs::empty, hash_uintE> H; //Set of elements such that deg(x) >= |H|
  size_t hindex; //|H|

  sparse_table<uintE, pbbs::empty, hash_uintE> P;

  pbbslib::dyn_arr<uintE> B; //Set of elements such that deg(x) == |H|

  sequence<pbbslib::dyn_arr<uintE>*> lowDegC; //Map of elements not in H sorted by degree up to a threshold
  sparse_table<uintE, pbbslib::dyn_arr<uintE>*, hash_uintE> highDegC; // Map of elements not in H sorted by degree greater than the threshold
  uintE threshold;

  // Makes it easier for sparse_table.find() which needs a default value
  pbbslib::dyn_arr<uintE>* empty;



  //a - threshold for low and high degree vertices
  HSet(uintE a, symmetric_graph<symmetric_vertex, pbbs::empty> _G) {

    *empty = pbbslib::dyn_arr<uintE>(0);

    *G = _G;
    threshold = a;


    H = make_sparse_table<uintE, pbbs::empty, hash_uintE> //|H| has to be between m/n and sqrt(2m)
      (std::max(G->m/G->n, (size_t) sqrt(2 * G->m)), std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());
    hindex = 0;

    P = make_sparse_table<uintE, pbbs::empty, hash_uintE>
      (std::max(G->m/G->n, (size_t) sqrt(2 * G->m)), std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());

    B = pbbslib::dyn_arr<uintE>(G->n);

    // C = make_sparse_table<uintE, sparse_table<uintE, pbbs::empty, hash_uintE>, hash_uintE>
      // (n, std::make_tuple(UINT_E_MAX, empty), hash_uintE());

    lowDegC = sequence<pbbslib::dyn_arr<uintE>*>(a); //possible degrees range from 0 ... a - 1
    highDegC = make_sparse_table<uintE, pbbslib::dyn_arr<uintE>*, hash_uintE>
      (G->n, std::make_tuple(UINT_E_MAX, new pbbslib::dyn_arr<uintE>(0)), hash_uintE());

  }

  uintE insert(uintE v) {

    //ADD TO GRAPH!!
    auto deg = G->get_vertex(v).getOutDegree();
    addToC(v);

    if(deg > hindex) {

      H.insert(std::make_tuple(v, pbbs::empty()));
      hindex++;

      if(B.size != 0) {
        auto y = B.A[0];
        B.erase(y);
        H.erase(y);
        hindex--;

        addToC(y);
      }

      else {
        // If C has entry for |H|, set that to B. Otherwise B is empty
        if (containedInC(hindex)) {
          
          if (hindex > threshold) {
            B = *highDegC.find(hindex, empty);
            highDegC.erase(hindex);
          }
          else {
            B = *lowDegC[hindex];
            lowDegC[hindex]->clear();
          }
        }

        else {
          B.clear();
        }

      }
    }

    return hindex;

  }

   
  uintE erase(uintE v) {

    // Remove from Graph!!

    auto deg = G->get_vertex(v).getOutDegree();
    if (B.indexOf(v) != -1) {
      B.erase(v);
    }

    else {
      removeFromC(v);
    }
    
    if (H.contains(v)) {

      auto h = hindex; //h is |H| before removing v
      H.erase(v);
      hindex--;

      if (P.contains(v)) {
        P.erase(v);
      }

      if (h > threshold) {

       if (highDegC.find(h, empty)->size != 0) { //Has entry for C[h]
          auto y = highDegC.find(h, empty)->A[0];
          removeFromC(y);
          B.add(y);
          H.insert(std::make_tuple(y, pbbs::empty()));
        }
        else {
          //highDegC.insert(std::make_tuple(h, *B));
          //B.clear(); //Might change the value of the entry in C 

          highDegC.insert(std::make_tuple(h, new pbbslib::dyn_arr(B))); //Creates pointer to copy of B
          B.clear();
        }
      }
      else {
        if (lowDegC[h] != nullptr || lowDegC[h]->size != 0) { //Has entry for C[h]
          auto y = lowDegC[h]->A[0];
          removeFromC(y);
          B.add(y);
          H.insert(std::make_tuple(y, pbbs::empty()));
        }
        else {
          lowDegC[h] = new pbbslib::dyn_arr(B); //Creates pointer to copy of B
          B.clear();
        }

      }
    }

    return hindex;
  }
	
  uintE change(uintE v) {

    erase(v);
    insert(v);

    //Why only add elements to P in an update?
    //What is the point of maintaining the dictionary parallel to C? Can't we just add during the update?
    if (G->get_vertex(v).getOutDegree() >= 2 * hindex) {
      P.insert(std::make_tuple(v, pbbs::empty()));
    }

    return hindex;
  }

  void addToC(uintE v) {

    uintE deg = G->get_vertex(v).getOutDegree();
    auto vertices = highDegC.find(deg, empty); //Vertices of degree "deg"

    // Add to highDegC
    if (deg > threshold) {
	    
      if (vertices->capacity == 0) { //If entry for deg is empty, create new sparse_table entry
          highDegC.insert(std::make_tuple(deg, new pbbslib::dyn_arr<uintE>(G->n)));
      }
      highDegC.find(deg, empty)->add(v);
    }
    // Add to lowDegC
    else {
      if (lowDegC[deg] == nullptr) {
        lowDegC[deg] = new pbbslib::dyn_arr<uintE>(G->n);
      }
      lowDegC[deg]->add(v);
    }

  }

  void removeFromC(uintE v) { // Assumes v is in C
    
    uintE deg = G->get_vertex(v).getOutDegree();
    
    //Remove from highDegC
    if (deg > threshold) {
      auto vertices = highDegC.find(deg, empty);
      if (vertices->size == 1) { // Deletes entire entry if it is empty (or will be)
        highDegC.erase(deg);
      }
      else {
        vertices->erase(v);
      }
    }
    //Remove from lowDegC
    else {

      lowDegC[deg]->erase(v);
    }
  }

  bool containedInC(uintE deg) {
    if (deg > threshold) {

      return highDegC.contains(deg);
    }
    else {
      if (lowDegC[deg] == nullptr) {
        return false;
      }
      return lowDegC[deg]->size != 0;
    }
    
  }

};


struct graph {

  public:
    HSet* h;
    // TODO: don't use constants like these; should be malloc-ed
    // using sequence would probably be easier than using malloc directly
    // Also, edges should come from the passed-in Graph template object
    std::vector<int> edges[1000];
    std::vector<int> hedges[100][100];

    int triangles = 0;

    graph(symmetric_graph<symmetric_vertex, pbbs::empty> _G) {
      h = new HSet(100, _G);
    }

    void insertV(uintE v) {
      h->insert(v);
    }

    int connect(uintE v, uintE u) {
      // Add u and v to H
      //h.inc(v);
      //h.inc(u);

      triangles += hedges[u][v].size();
      for(int i = 0;i < edges[u].size();i++) {
        // Add new wedge v -- u -- i
        hedges[v][edges[u][i]].push_back(u);
        hedges[edges[u][i]][v].push_back(u);
      }
      for(int i = 0;i < edges[v].size();i++) {
        // Add new wedge u -- v -- i
        hedges[u][edges[v][i]].push_back(v);
        hedges[edges[v][i]][u].push_back(v);
      }
      edges[v].push_back(u);
      edges[u].push_back(v);

      return triangles;
    }
};
