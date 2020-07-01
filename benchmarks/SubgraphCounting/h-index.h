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
#include "pbbslib/sequence.h"

template <class Graph>
struct HSet {

  Graph* G; //Graph

  sparse_table<uintE, pbbs::empty, hash_uintE> H; //Set of elements such that deg(x) >= |H|
  size_t hindex; //|H|

  sparse_table<uintE, pbbs::empty, hash_uintE> P;

  uintE B; //Number of elements such that deg(x) == |H|

  pbbs::sequence<pbbs::sequence<uintE>*> C; //Map of elements not in H fed by degree up to a threshold

  //Doesn't have to be symmetric vertex
  //Specify templatized graph
  //a - threshold for low and high degree vertices
  HSet(Graph& _G) {

    G = &_G;

    
    H = make_sparse_table<uintE, pbbs::empty, hash_uintE> //|H| has to be between m/n and sqrt(2m)
      (std::max(G->m/G->n, (size_t) sqrt(2 * G->m)), std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());
    hindex = 0;

    P = make_sparse_table<uintE, pbbs::empty, hash_uintE>
      (std::max(G->m/G->n, (size_t) sqrt(2 * G->m)), std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());

    B = 0;

    C = sequence<pbbs::sequence<uintE>*>(G->n + 10); //possible degrees range from 0 ... a
    par_for(0, C.size(), [&] (size_t i) { C[i] = new sequence<uintE>(0); });
  }

  //--------------------------INSERT--------------------------//
  //batch is a sequence of vertex ids (that can be used to get from graph)
  uintE insert(sequence<uintE> batch, bool sorted = false) {

    pbbs::sequence<uintE> sortedBatch;
    if (!sorted) {
      sortedBatch = integer_sort(batch, [&] (uintE v) { return G->get_vertex(v).getOutDegree(); }).rslice();
    }
    else {
      sortedBatch = batch;
    } 
    batch.clear();
    //Stores degrees and sorts in descending order
    sequence<uintE> deg = pbbs::sequence<uintE>(sortedBatch.size());
    par_for(0, sortedBatch.size(), [&] (size_t i) {
      deg[i] = deg[i] = G->get_vertex(sortedBatch[i]).getOutDegree();
    });
    
    //--------------------------Compute H-index--------------------------//
    pbbs::sequence<uintE> sum = pbbs::sequence<uintE>(sortedBatch.size());
    
    par_for(0, sortedBatch.size(), [&] (size_t i) {
      sum[i] = C[hindex + i + 1]->size();
    });
    //Computes Prefix Sum (note scan_add_inplace is exclusive)
    pbbslib::scan_add_inplace(sum);
    par_for(0, sum.size(), [&] (size_t i) {
      sum[i] += C[hindex + i + 1]->size() + B; //Compute inclusive sunm
    });
    
    //Creates the binary array described as M
    pbbs::sequence<bool> M = pbbs::sequence<bool>(sum.size() + 1);
    if (B < deg.size()) {
      M[0] = (deg[B] > hindex);
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
    uintE greaterThanH = hindex + (C[hindex]->size() - B); //Number of vertices greater than hindex
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
    B = C[hindex]->size() - ((aboveNewH + added) - hindex); //Number of vertices not in B is equal to the number of vertices greater than hindex minus hindex

    //Slight variation to hindex paper - C will store all vertices (instead of the vertices not in B)


    deg.clear();
    sum.clear();   

    return hindex;
  }
  
  //--------------------------ERASE--------------------------//   
  uintE erase(sequence<uintE> batch, bool sorted = false) {

    pbbs::sequence<uintE> sortedBatch;
    if (!sorted) {
      sortedBatch = integer_sort(batch, [&] (uintE v) { return G->get_vertex(v).getOutDegree(); }).rslice();
    }
    else {
      sortedBatch = batch;
    } 
    batch.clear();
    //Stores degrees and sorts in ascending order
    sequence<uintE> deg = pbbs::sequence<uintE>(sortedBatch.size());
    par_for(0, sortedBatch.size(), [&] (size_t i) {
      deg[i] = G->get_vertex(sortedBatch[i]).getOutDegree();
    });

    //--------------------------Compute H-index--------------------------//
    pbbs::sequence<uintE> sum = pbbs::sequence<uintE>(std::min(sortedBatch.size(), hindex) + 1).rslice();
    
    par_for(0, sum.size(), [&] (size_t i) {
      sum[i] = C[hindex - i]->size();
    });

    //Computes Prefix Sum (note scan_add_inplace is exclusive)
    pbbslib::scan_add_inplace(sum);
    par_for(0, sum.size(), [&] (size_t i) {
      sum[i] += C[hindex - i]->size(); //Compute inclusive sunm
    });

    pbbs::sequence<bool> M = pbbs::sequence<bool>(sum.size());

    //Creates the binary array described as M
    if (C[hindex]->size() - B < deg.size()) {
      M[0] = (deg[C[hindex]->size() - B] >= hindex);
    }
    else {
      M[0] = false;
    }
    par_for(1, sum.size(), [&] (size_t i) {
      if (sum[i] - i < deg.size()) {
        M[i] = (deg[sum[i] + i] >= hindex - i);
      }
      else {
        M[i] = false;
      }
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
    
    //Updating B
    uintE greaterThanH = hindex + (C[hindex]->size() - B); //Number of vertices greater than hindex
    pbbs::sequence<uintE> sizesBetweenH = pbbs::sequence<uintE>(hindexDecrease); //Sizes of each entry in C that is between old hindex
    par_for(0, hindexDecrease, [&] (size_t i) {
      sizesBetweenH[i] = C[hindex - i - 1]->size();
    });
    uintE betweenH = pbbslib::reduce_add(sizesBetweenH); //Number of vertices between old hindex and (new hindex - 1)
    uintE aboveNewH = greaterThanH + betweenH;

    uintE removed = pbbs::binary_search(deg, [&] (uintE x) { //Number of vertices that are above hindex
      return x >= hindex - hindexDecrease;
    });
    
    hindex -= hindexDecrease;
    //Update C
    removeFromC(sortedBatch, deg);
    B = C[hindex]->size() - ((aboveNewH - removed) - hindex); //Number of vertices not in B is equal to the number of vertices greater than hindex minus hindex

    return hindex;
  }
  

  //--------------------------ADD TO C--------------------------//
  void addToC(pbbs::sequence<uintE> batch, pbbs::sequence<uintE> deg) {
    
    //Subraction
    sequence<bool> difference = pbbs::sequence<bool>(batch.size() - 1);
    
    par_for(0, deg.size() - 1, [&] (size_t i) {
      difference[i] = (deg[i] != deg[i + 1]);  
    });

    sequence<size_t> idx = pbbs::sequence<size_t>(difference.size() + 1);

    par_for(0, difference.size() + 1, [&] (size_t i) {
      idx[i] = i;
    });

    auto f = [&] (size_t i) { return difference[i] || i == batch.size() - 1; };
    auto indices = pbbs::filter(idx, f);
    idx.clear();
    difference.clear();
    //Uses indices sequence to know the index of last entry for each clustered deg
    par_for(0, indices.size(), [&] (size_t i) {
      size_t start = (i == 0 ? 0 : indices[i - 1] + 1);
      size_t end = indices[i] + 1;

      pbbs::sequence<uintE> extra = batch.slice(start, end);
      
      C[deg[indices[i]]]->append(extra);
    });
    indices.clear();
  }
  
  //--------------------------REMOVE FROM C--------------------------//
  void removeFromC(pbbs::sequence<uintE> batch, pbbs::sequence<uintE> deg) {
    
    //Subraction
    sequence<bool> difference = pbbs::sequence<bool>(batch.size() - 1);
    
    par_for(0, deg.size() - 1, [&] (size_t i) {
      difference[i] = (deg[i] != deg[i + 1]);  
    });

    sequence<size_t> indices = pbbs::sequence<size_t>(difference.size() + 1);

    par_for(0, difference.size() + 1, [&] (size_t i) {
      indices[i] = i;
    });
    // TODO: Why do you need a separate difference array to make this function f?
    // Can't you just check deg directly?
    auto f = [&] (size_t i) { return difference[i] || i == batch.size() - 1; };
    indices = pbbs::filter(indices, f);
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
      auto temp = pbbs::filter(*C[deg[indices[i]]], remove);
      C[deg[indices[i]]] = &temp;
    });
  }
  
  //batch[i] = vertex id, degs[i] = vertex degree
  uintE change(sequence<uintE> batch, sequence<uintE> degs) {
    pbbs::sequence<uintE> sortedBatch = integer_sort(batch, [&] (uintE v) { return G->get_vertex(v).getOutDegree(); }).rslice();
    erase(sortedBatch, true);
    //Change in graph using batch and degs
    insert(sortedBatch, true);

    return hindex;
  }

  //Returns tuple of contains and the sparse table of updated H
  //Instead of iterating through everything in insertion or deletion, we could just do it inside triangle counting
  //When we test for elements inside H, for the first vertex we can iterate through C looking for a vertex while adding everything to H
  //For the rest of the vertices, we can just use H.contains(target)
  std::tuple<bool, sparse_table<uintE, pbbs::empty, hash_uintE>> containedInH(uintE target) {
    size_t num = 0;
    uintE deg = hindex;
    bool flag = false;
    while (num < hindex) {

      if (deg == hindex) {
        for (size_t i = 0; i < B; i++) {

          uintE v = (*C[deg])[i];
          H.insert(std::make_tuple(v, pbbs::empty()));
          num++;

          if (v == target) {
            flag = true;
          }
        }
        deg++;
      }

      else {

        for (size_t i = 0; i < C[deg]->size(); i++) {

          uintE v = (*C[deg])[i];
          H.insert(std::make_tuple(v, pbbs::empty()));
          num++;

          if (v == target) {
            flag = true;
          }
        }
        deg++;
      }
    }

    return std::make_tuple(flag, H);
  }

};

/*
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
      h = new HSet(_G);
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
*/
