//- Return an iterator of H
//  - Return a lambda
//  - An array as a function

#include <iostream>
#include <map>
#include <set>
#include <limits>
#include <vector>
#include <math.h>
#include "ligra/ligra.h"
#include "dynamic_symmetric_graph.h"
#include "ligra/pbbslib/sparse_table.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "pbbslib/sequence.h"

template <class Graph>
struct HSet {

  //Graph* G; //Graph
  dynamic_symmetric_graph<dynamic_symmetric_vertex, pbbs::empty>* G;

  sparse_table<uintE, pbbs::empty, hash_uintE> H; //Set of elements such that deg(x) >= |H|
  size_t hindex; //|H|

  //sparse_table<uintE, pbbs::empty, hash_uintE> P;

  //uintE bSize; //Number of elements in B
  pbbs::sequence<uintE> B;

  pbbs::sequence<pbbs::sequence<uintE>*> C; //Map of elements not in H fed by degree up to a threshold


  //Constructor
  //HSet(Graph& _G) { //Don't pass by reference
  HSet(dynamic_symmetric_graph<dynamic_symmetric_vertex, pbbs::empty>* _G) {
    G = _G;
    
    H = make_sparse_table<uintE, pbbs::empty, hash_uintE> //|H| has to be between m/n and sqrt(2m)
      (std::max(G->m/G->n, (size_t) sqrt(2 * G->m)), std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());
    hindex = 0;

    B = pbbs::sequence<uintE>(0);

    C = pbbs::sequence<pbbs::sequence<uintE>*>(G->n);
    par_for(0, C.size(), [&] (size_t i) { C[i] = new sequence<uintE>(0); });

    //Add all starting graph elements into HSet
    auto start = pbbs::sequence<uintE>(G->n);
    par_for(0, start.size(), [&] (size_t i) {
      start[i] = i;
    });
    insert(start);
  }

  //--------------------------INSERT--------------------------//
  //batch is a sequence of vertex ids (that can be used to get from graph)
  //sorted = true if batch is already sorted in descending order
  uintE insert(sequence<uintE> batch, bool sorted = false) {

    pbbs::sequence<uintE> sortedBatch;
    if (!sorted) {
      sortedBatch = integer_sort(batch, [&] (uintE v) { return G->get_vertex(v).degree; }).rslice();
    }
    else {
      sortedBatch = batch;
    }
    batch.clear();
    
    //Stores degrees and sorts in descending order
    sequence<uintE> deg = pbbs::sequence<uintE>(sortedBatch.size());
    par_for(0, sortedBatch.size(), [&] (size_t i) {
      deg[i] = deg[i] = G->get_vertex(sortedBatch[i]).degree;
    });

    //Adds additional space to C if it is too small
    if (sortedBatch.size() >= C.size()) { //deg[0] is the largest degree
      
      size_t cOldSize = C.size();

      //Double it or increase to the size of batch (no errors when slicing)
      C.resize(std::max(C.size(), sortedBatch.size() - C.size()), nullptr);
      par_for(cOldSize, C.size(), [&] (size_t i) { 
        C[i] = new sequence<uintE>(0); 
      });
    }
    
    
    //--------------------------Compute H-index--------------------------//
    pbbs::sequence<uintE> sum = pbbs::sequence<uintE>(sortedBatch.size());

    par_for(0, sortedBatch.size(), [&] (size_t i) {
      sum[i] = C[hindex + i + 1]->size();
    });

    //Computes Prefix Sum (note scan_add_inplace is exclusive)
    pbbslib::scan_add_inplace(sum);
    par_for(0, sum.size(), [&] (size_t i) {
      sum[i] += C[hindex + i + 1]->size() + B.size(); //Compute inclusive sunm
    });

    //Creates the binary array described as M
    pbbs::sequence<bool> M = pbbs::sequence<bool>(sum.size() + 1);
    if (B.size() < deg.size()) {
      M[0] = (deg[B.size()] > hindex);
    }
    else {
      M[0] = false;
    }
    par_for(1, sum.size() + 1, [&] (size_t i) {
      if (sum[i - 1] + i < deg.size()) {
        M[i] = (deg[sum[i - 1] + i] > hindex + i);
      }
      else {
        M[i] = false;
      }
    });

    //Find index of first 0
    pbbs::sequence<size_t> indexM = pbbs::sequence<size_t>(M.size());
    par_for(0, M.size(), [&] (size_t i) {
      indexM[i] = !M[i] ? i : UINT_E_MAX;
    });

    uintE hindexIncrease = pbbslib::reduce_min(indexM);
    M.clear();
    indexM.clear();

    if (hindexIncrease == UINT_E_MAX) {
      hindexIncrease = sortedBatch.size();
    }
    //Update B
    uintE greaterThanH = hindex + (C[hindex]->size() - B.size()); //Number of vertices greater than hindex
    pbbs::sequence<uintE> sizesBetweenH = pbbs::sequence<uintE>(hindexIncrease); //Sizes of each entry in C that is between old hindex
    par_for(0, hindexIncrease, [&] (size_t i) {
      sizesBetweenH[i] = C[hindex + i]->size();
    });
    uintE betweenH = pbbslib::reduce_add(sizesBetweenH); //Number of vertices between old hindex and (new hindex - 1)
    uintE aboveNewH = greaterThanH - betweenH;

    uintE added = pbbs::binary_search(deg, [&] (uintE x) { //Number of vertices that are above hindex
      return x >= hindex + hindexIncrease;
    });

    addToC(sortedBatch, deg);
    hindex += hindexIncrease;
    //Update C, hindexIncrease);
    uintE bSize = C[hindex]->size() - ((aboveNewH + added) - hindex); //Number of vertices not in B is equal to the number of vertices greater than hindex minus hindex
    B = C[hindex]->slice(0, bSize);
    //Slight variation to hindex paper - C will store all vertices (instead of the vertices not in B)

    deg.clear();
    sum.clear();

    return hindex;
  }
  
  //--------------------------ERASE--------------------------//   
  //batch is a sequence of vertex ids (that can be used to get from graph)
  //sorted = true if batch is already sorted in descending order
  uintE erase(sequence<uintE> batch, bool sorted = false) {


    pbbs::sequence<uintE> sortedBatch;
    if (!sorted) {
      sortedBatch = integer_sort(batch, [&] (uintE v) { return G->get_vertex(v).degree; }).rslice();
    }
    else {
      sortedBatch = batch;
    } 
    batch.clear();
    //Stores degrees and sorts in ascending order
    sequence<uintE> deg = pbbs::sequence<uintE>(sortedBatch.size());
    par_for(0, sortedBatch.size(), [&] (size_t i) {
      deg[i] = G->get_vertex(sortedBatch[i]).degree;
    });


    //--------------------------Updating Variables--------------------------//

    //Set up things to find new B
    uintE greaterThanH = hindex + (C[hindex]->size() - B.size()); //Number of vertices greater than hindex

    pbbs::sequence<uintE> sizesBetweenH = pbbs::sequence<uintE>(std::min(sortedBatch.size(), hindex)); //Sizes of each entry in C that is between old hindex
    par_for(0, sizesBetweenH.size(), [&] (size_t i) {
      sizesBetweenH[i] = C[hindex - i - 1]->size();
    });

    //Number of elements in H that were removed
    uintE removed = pbbs::binary_search(deg, [&] (uintE x) { //Number of vertices that are above hindex
      return x >= hindex;
    });

    uintE stillAboveH = greaterThanH - removed;


    //Update C
    removeFromC(sortedBatch, deg);
    //Update Current B
    auto remove = [&] (uintE element) {
      for (size_t i = 0; i < sortedBatch.size(); i++) {
        if (element == sortedBatch[i]) return false;
      }
      return true;
    };

    //Copies a lot which might be bad
    //Could also modify filter to return pointer
    auto filtered = pbbs::filter(B, remove);
    B.clear();
    B = pbbs::sequence<uintE>(filtered);
    filtered.clear();


    //--------------------------Compute H-index--------------------------//
    pbbs::sequence<uintE> sum = pbbs::sequence<uintE>(std::min(sortedBatch.size(), hindex)).rslice();
    
    par_for(0, sum.size(), [&] (size_t i) {
      sum[i] = C[hindex - i - 1]->size();
    });

    //Computes Prefix Sum (note scan_add_inplace is exclusive)
    pbbslib::scan_add_inplace(sum);
    par_for(0, sum.size(), [&] (size_t i) {
      sum[i] += C[hindex - i]->size(); //Compute inclusive sunm
    });

    pbbs::sequence<bool> M = pbbs::sequence<bool>(sum.size() + 1);

    M[0] = (stillAboveH < hindex);
    par_for(1, sum.size(), [&] (size_t i) {
      M[i] = (sum[i - 1] + stillAboveH < hindex - i);
    });
    sum.clear();

    //Find index of first 0
    pbbs::sequence<size_t> indexM = pbbs::sequence<size_t>(M.size());
    par_for(0, M.size(), [&] (size_t i) {
      indexM[i] = !M[i] ? i : UINT_E_MAX;
    });

    uintE hindexDecrease = pbbslib::reduce_min(indexM);
    M.clear();
    indexM.clear();

    if (hindexDecrease == UINT_E_MAX) {
      hindexDecrease = std::min(sortedBatch.size(), hindex);
    }

    //Slight variation to hindex paper - C will store all vertices (instead of the vertices not in B)
    pbbs::sequence<uintE> sizesBetween = sizesBetweenH.slice(0, hindexDecrease);
    sizesBetweenH.clear();

    uintE betweenH = pbbslib::reduce_add(sizesBetween); //Number of vertices between old hindex and (new hindex - 1)
    sizesBetween.clear();

    uintE aboveNewH = greaterThanH + betweenH;

    uintE removedFromNewH = pbbs::binary_search(deg, [&] (uintE x) { //Number of vertices that are above hindex
      return x >= hindex - hindexDecrease;
    });

    hindex -= hindexDecrease;

    //Number of vertices not in B is equal to the number of vertices greater than hindex minus hindex
    uintE bSize = C[hindex]->size() - ((aboveNewH - removedFromNewH) - hindex);
    B = C[hindex]->slice(0, bSize);

    return hindex;
  }
  

  //--------------------------ADD TO C--------------------------//
  void addToC(pbbs::sequence<uintE> batch, pbbs::sequence<uintE> deg) {

    /*    
    //Subraction
    sequence<bool> difference = pbbs::sequence<bool>(batch.size() - 1);
    
    par_for(0, deg.size() - 1, [&] (size_t i) {
      difference[i] = (deg[i] != deg[i + 1]);  
    });
    */
    if (deg[0] >= C.size()) {
      size_t cOldSize = C.size();

      //Double it or increase it up to the largest degree
      C.resize(std::max(C.size(), deg[0] - C.size()), nullptr);
      par_for(cOldSize, C.size(), [&] (size_t i) {
        C[i] = new sequence<uintE>(0); 
      });
    }

    sequence<size_t> idx = pbbs::sequence<size_t>(batch.size());

    par_for(0, batch.size(), [&] (size_t i) {
      idx[i] = i;
    });

    auto f = [&] (size_t i) { return i == batch.size() - 1 || deg[i] != deg[i + 1]; };
    auto indices = pbbs::filter(idx, f);
    idx.clear();
    //difference.clear();
    //Uses indices sequence to know the index of last entry for each clustered deg
    par_for(0, indices.size(), [&] (size_t i) {
      size_t start = (i == 0 ? 0 : indices[i - 1] + 1);
      size_t end = indices[i] + 1;

      pbbs::sequence<uintE> extra = batch.slice(start, end);
      
      C[deg[indices[i]]]->append(extra);
      extra.clear();
    });
    indices.clear();
  }
  
  //--------------------------REMOVE FROM C--------------------------//
  void removeFromC(pbbs::sequence<uintE> batch, pbbs::sequence<uintE> deg) {
    
    /*
    //Subraction
    sequence<bool> difference = pbbs::sequence<bool>(batch.size() - 1);
    
    par_for(0, deg.size() - 1, [&] (size_t i) {
      difference[i] = (deg[i] != deg[i + 1]);  
    });
    */

    sequence<size_t> idx = pbbs::sequence<size_t>(batch.size());

    par_for(0, batch.size(), [&] (size_t i) {
      idx[i] = i;
    });

    auto f = [&] (size_t i) { return i == batch.size() - 1 || deg[i] != deg[i + 1]; };
    auto indices = pbbs::filter(idx, f);
    //Uses indices sequence to know the index of last entry for each clustered deg
    par_for(0, indices.size(), [&] (size_t i) {
      size_t start = (i == 0 ? 0 : indices[i - 1] + 1);
      size_t end = indices[i] + 1;

      pbbs::sequence<uintE> extra = batch.slice(start, end);
      
      auto remove = [&] (uintE element) {
        for (size_t i = 0; i < extra.size(); i++) {
          if (element == extra[i]) return false;
        }
        return true;
      };

      //Copies a lot which might be bad
      //Could also modify filter to return pointer
      auto filtered = pbbs::filter(*C[deg[indices[i]]], remove);
      C[deg[indices[i]]]->clear();
      C[deg[indices[i]]] = new pbbs::sequence<uintE>(filtered);
      filtered.clear();
    });
    indices.clear();
  }
  
  //Modifies the degree of existing vertices
  //batch[i] = vertex id, degs[i] = vertex degree
  //--------------------------CHANGE--------------------------//
  uintE change(sequence<uintE> batch, sequence<uintE> degs) {
    pbbs::sequence<uintE> sortedBatch = integer_sort(batch, [&] (uintE v) { return G->get_vertex(v).degree; }).rslice();
    erase(sortedBatch, true);
    //Change in graph using batch and degs
    insert(sortedBatch, true);

    return hindex;
  }

  //Contains can just check the degree of a vertex and see if it is greater than |H|
  //If vertex's degree is the hindex, we check if it is contained in B
  bool contains(uintE target) {
    uintE deg = G->get_vertex(target).degree;

    if (deg > hindex) return true;
    else if (deg == hindex) {
      for (size_t i = 0; i < B.size(); i++) {
        if (B[i] == target) return true;
      }
    }
    return false;
  }

};
