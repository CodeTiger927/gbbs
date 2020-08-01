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


size_t INIT_DEG_SIZE = 4;

template <class Graph>
struct HSet {

  //Graph* G; //Graph
  dynamic_symmetric_graph<dynamic_symmetric_vertex, pbbs::empty>* G;

  size_t hindex; //|H|

  //sparse_table<uintE, pbbs::empty, hash_uintE> P;

  uintE bSize; //Number of elements in B
  //Can access B by finding first bSize elements inside C

  //Switch to dyn_arr
  //Store loose upper bound

  pbbslib::dyn_arr<pbbslib::dyn_arr<uintE>*> C; //Map of elements not in H fed by degree up to a threshold


  //Constructor
  //HSet(Graph& _G) { //Don't pass by reference
  HSet(dynamic_symmetric_graph<dynamic_symmetric_vertex, pbbs::empty>* _G) {
    G = _G;
    
    hindex = 0;

    bSize = 0;

    C = pbbslib::dyn_arr<pbbslib::dyn_arr<uintE>*>(100); //FIX LATER
    par_for(0, C.capacity, [&] (size_t i) { C.A[i] = nullptr; });

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
    
    //--------------------------Compute H-index--------------------------//
    pbbs::sequence<uintE> sum = pbbs::sequence<uintE>(sortedBatch.size());

    par_for(0, sortedBatch.size(), [&] (size_t i) {
      if (hindex + i + 1 < C.capacity && C.A[hindex + i + 1] != nullptr) {
        sum[i] = C.A[hindex + i + 1]->size;
      }
      else {
        sum[i] = 0;
      }
    });

    //Computes Prefix Sum (note scan_add_inplace is exclusive)
    pbbslib::scan_add_inplace(sum);
    par_for(0, sum.size(), [&] (size_t i) {
      if (hindex + i + 1 < C.capacity && C.A[hindex + i + 1] != nullptr) {
        sum[i] += C.A[hindex + i + 1]->size + bSize; //Compute inclusive sum
      }
      else {
        sum[i] += bSize;
      }
    });

    //Creates the binary array described as M
    pbbs::sequence<bool> M = pbbs::sequence<bool>(sum.size() + 1);
    if (bSize < deg.size()) {
      M[0] = (deg[bSize] > hindex);
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

    uintE greaterThanH;
    if (C.A[hindex] != nullptr) {
      greaterThanH = hindex + (C.A[hindex]->size - bSize); //Number of vertices greater than hindex
    }
    else {
      greaterThanH = hindex + (0 - bSize);
    }
    
    pbbs::sequence<uintE> sizesBetweenH = pbbs::sequence<uintE>(hindexIncrease); //Sizes of each entry in C that is between old hindex
    par_for(0, hindexIncrease, [&] (size_t i) {
      if (hindex + i < C.capacity && C.A[hindex + i] != nullptr) {
        sizesBetweenH[i] = C.A[hindex + i]->size; 
      }
      else {
        sizesBetweenH[i] = 0;
      }
    });
    uintE betweenH = pbbslib::reduce_add(sizesBetweenH); //Number of vertices between old hindex and (new hindex - 1)
    uintE aboveNewH = greaterThanH - betweenH;

    uintE added = pbbs::binary_search(deg, [&] (uintE x) { //Number of vertices that are above hindex
      return x >= hindex + hindexIncrease;
    });

    addToC(sortedBatch, deg);
    hindex += hindexIncrease;
    //Update C, hindexIncrease);
    if (C.A[hindex] != nullptr) {
      //Number of vertices not in B is equal to the number of vertices greater than hindex minus hindex
      bSize = C.A[hindex]->size - ((aboveNewH + added) - hindex); 
    }
    else {
      bSize = -((aboveNewH + added) - hindex);
    }
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
    uintE greaterThanH;
    if (C.A[hindex] != nullptr) {
      greaterThanH = hindex + (C.A[hindex]->size - bSize); //Number of vertices greater than hindex
    }
    else {
      greaterThanH = hindex + (0 - bSize);
    }

    pbbs::sequence<uintE> sizesBetweenH = pbbs::sequence<uintE>(std::min(sortedBatch.size(), hindex)); //Sizes of each entry in C that is between old hindex
    par_for(0, sizesBetweenH.size(), [&] (size_t i) {
      if (C.A[hindex - i - 1] != nullptr) {
        sizesBetweenH[i] = C.A[hindex - i - 1]->size;
      }
      else {
        sizesBetweenH[i] = 0;
      }
    });

    //Number of elements in H that were removed
    uintE removed = pbbs::binary_search(deg, [&] (uintE x) { //Number of vertices that are above hindex
      return x >= hindex;
    });

    uintE stillAboveH = greaterThanH - removed;

    //Update C
    removeFromC(sortedBatch, deg);


    //--------------------------Compute H-index--------------------------//
    pbbs::sequence<uintE> sum = pbbs::sequence<uintE>(std::min(sortedBatch.size(), hindex)).rslice();
    
    par_for(0, sum.size(), [&] (size_t i) {
      if (C.A[hindex - i - 1] != nullptr) {
        sum[i] = C.A[hindex - i - 1]->size;
      }
      else {
        sum[i] = 0;
      }
    });

    //Computes Prefix Sum (note scan_add_inplace is exclusive)
    pbbslib::scan_add_inplace(sum);
    par_for(0, sum.size(), [&] (size_t i) {
      if (hindex - i - 1 >= 0 && C.A[hindex - i - 1] != nullptr) {
        sum[i] += C.A[hindex - i - 1]->size; //Compute inclusive sum
      }
      else {
        sum[i] += 0;
      }
    });

    pbbs::sequence<bool> M = pbbs::sequence<bool>(sum.size() + 1);

    M[0] = (stillAboveH < hindex);
    par_for(1, M.size(), [&] (size_t i) {
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
    if (C.A[hindex] != nullptr) {
      bSize = C.A[hindex]->size - ((aboveNewH - removedFromNewH) - hindex);
    }
    else {
      bSize = -((aboveNewH - removedFromNewH) - hindex);
    }

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
      
      if (C.A[deg[indices[i]]] == nullptr) {
        C.A[deg[indices[i]]] = new pbbslib::dyn_arr<uintE>(INIT_DEG_SIZE);
      }
      C.A[deg[indices[i]]]->add(extra);
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

      if (deg[indices[i]] != hindex) {
        auto filtered = pbbs::new_array_no_init<uintE>(C.A[deg[indices[i]]]->capacity);
        size_t removed = pbbslib::filter_seq(C.A[deg[indices[i]]]->A, filtered, C.A[deg[indices[i]]]->size, remove);
        C.A[deg[indices[i]]]->del();
        C.A[deg[indices[i]]]->A = filtered;
        C.A[deg[indices[i]]]->size = removed;
      }
      else {
        auto filtered = pbbs::new_array_no_init<uintE>(C.A[deg[indices[i]]]->capacity);

        //Temporarily puts things into filtered (will be replaced in next line)
        uintE removedFromB = pbbslib::filter_seq(C.A[deg[indices[i]]]->A, filtered, bSize, remove); 
        size_t removed = pbbslib::filter_seq(&(C.A[deg[indices[i]]]->A)[bSize], &(filtered)[removedFromB], C.A[deg[indices[i]]]->size - bSize, remove);
        C.A[deg[indices[i]]]->del();
        C.A[deg[indices[i]]]->A = filtered;
        C.A[deg[indices[i]]]->size = removed;

        bSize = removedFromB;
      }

      if (C.A[deg[indices[i]]]->size == 0) {
        C.A[deg[indices[i]]] = nullptr;
      }

      extra.clear();
    });

    indices.clear();
  }

  uintE insertVertices(sequence<uintE> vertices) {
    G->batchAddVertices(vertices);
    return insert(vertices);
  }

  uintE eraseVertices(sequence<uintE> vertices) {
    G->batchRemoveVertices(vertices);
    return erase(vertices);
  }
  
  //Insert edges once e.g. u--v inserts v--u as well
  uintE insertEdges(sequence<std::pair<uintE, uintE>> edges) {
    sequence<uintE> vertices = sequence<uintE>(2 * edges.size());
    par_for(0, edges.size(), [&] (size_t i) {
     vertices[2 * i] = edges[i].first;
     vertices[(2 * i) + 1] = edges[i].second; 
    });

    erase(vertices);
    G->batchAddEdges(edges);
    return insert(vertices);
  }

  uintE eraseEdges(sequence<std::pair<uintE, uintE>> edges) {
    sequence<uintE> vertices = sequence<uintE>(2 * edges.size());
    par_for(0, edges.size(), [&] (size_t i) {
      vertices[2 * i] = edges[i].first;
      vertices[(2 * i) + 1] = edges[i].second;
    });

    erase(vertices);
    G->batchRemoveEdges(edges);
    return insert(vertices);
  }

  //Contains can just check the degree of a vertex and see if it is greater than |H|
  //If vertex's degree is the hindex, we check if it is contained in B
  bool contains(uintE target) {
    uintE deg = G->get_vertex(target).degree;

    if (deg > hindex) return true;
    else if (deg == hindex) {
      for (size_t i = 0; i < bSize; i++) {
        if (C.A[hindex]->A[i] == target) return true;
      }
    }
    return false;
  }

};
