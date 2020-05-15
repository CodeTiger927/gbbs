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

// New Code for HSet
// For hashing sparse table

//Already defined
//struct hash_uintE {
  //inline size_t operator () (const uintE & a) {return pbbs::hash64_2(a);}
//};

//TODO: same notes as on h-index2


// TODO: we can talk about this more tonight, but I wasn't quite sure exactly what you
// were doing with the resizings and H set -- this is how I think the h-index code should go.
template<class Graph>
struct HSet2 {
  Graph* G = nullptr;
  pbbs::sequence<uintE*> C;
  pbbs::sequence<uintE> B;
  size_t h_index = 0;
  pbbs::sequence<uintE> old_degrees;

  HSet2(Graph* _G, pbbs::sequence<uintE> _old_degrees) : G(_G), old_degrees(_old_degrees) {}

  // Note: This batch is for vertices that move "up" -- their degrees increase only,
  // and H cannot decrease
  // We'll have another function for vertices whose degrees decrease only.
  uintE insert_additions(sequence<uintE>& batch) {
    // Note: G must already be updated with the new degrees of vertices in batch, but we should separately
    // keep an array with old degrees, that holds the degree state of the graph from before the update
    // First, remove the vertices in batch from C
      // Do this by looking up old degrees
    // Then, add the new degrees to C
      // Do this by making an array with each entry as deg(v) for v in batch, and do a sort
      // Then, subtract each element in the array from its previous element
      // Do a filter that compresses and returns the indices of the 1 elements
      // These indices tell you how many elements have each new degree, which we can use to figure out what should be increased
    // Filter from batch all elements with new degree <= |H|
    // From intersection of batch and B, remove these elements from both (degree must have increased)
    // Now, if |batch| <= |B|, remove |B| - |batch| elements from B
    // Otherwise, empty B and remove the lowest |B| elements from |batch|
       // The rest will increase our h index
       // To find out by how much, make an array with each entry as deg(v) - h_index for v in batch
       // Do a sort, then make a marking array where you mark any entry that's less than or equal to its index
       // Do a reduce for the minimum index that's marked -- this is how much our h-index is increased by
       // The entry in C for new h-index becomes our new B
  }

};





struct HSet {

  symmetric_graph<symmetric_vertex, pbbs::empty>* G; //Graph

  sparse_table<uintE, pbbs::empty, hash_uintE> H; //Set of elements such that deg(x) >= |H|
  size_t hindex; //|H|

  sparse_table<uintE, pbbs::empty, hash_uintE> P;

  pbbs::sequence<uintE> B; //Set of elements such that deg(x) == |H|

  pbbs::sequence<pbbs::sequence<uintE>*> lowDegC; //Map of elements not in H sorted by degree up to a threshold
  sparse_table<uintE, pbbs::sequence<uintE>*, hash_uintE> highDegC; // Map of elements not in H sorted by degree greater than the threshold
  uintE threshold;

  pbbs::sequence<uintE> lowDegD;
  sparse_table<uintE, uintE, hash_uintE> highDegD;


  //a - threshold for low and high degree vertices
  HSet(uintE a, symmetric_graph<symmetric_vertex, pbbs::empty> _G) {

    // TODO: Why can't this be nullptr?
    G = &_G;
    threshold = a;

    
    H = make_sparse_table<uintE, pbbs::empty, hash_uintE> //|H| has to be between m/n and sqrt(2m)
      (std::max(G->m/G->n, (size_t) sqrt(2 * G->m)), std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());
    hindex = 0;

    P = make_sparse_table<uintE, pbbs::empty, hash_uintE>
      (std::max(G->m/G->n, (size_t) sqrt(2 * G->m)), std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());

    B = pbbs::sequence<uintE>(0);

    // C = make_sparse_table<uintE, sparse_table<uintE, pbbs::empty, hash_uintE>, hash_uintE>
      // (n, std::make_tuple(UINT_E_MAX, empty), hash_uintE());

    lowDegC = sequence<pbbs::sequence<uintE>*>(a + 1); //possible degrees range from 0 ... a
    // TODO: You can reduce the search space of highDegC by taking G->n - threshold
    // It won't make a huge difference, but maybe it'll be better for the hash
    // function to distribute these values evenly, since we don't have any
    // values in ehre from 0 .. threshold
    highDegC = make_sparse_table<uintE, pbbs::sequence<uintE>*, hash_uintE>
      (G->n, std::make_tuple(UINT_E_MAX, nullptr), hash_uintE());

    par_for(0, lowDegC.size(), [&] (size_t i) { lowDegC[i] = new sequence<uintE>(0); });

    lowDegD = sequence<uintE>(a);
    par_for(0, lowDegD.size(), [&] (size_t i) { lowDegD[i] = 0; });
    highDegD = make_sparse_table<uintE, uintE, hash_uintE>
      (G->n, std::make_tuple(UINT_E_MAX, UINT_E_MAX), hash_uintE());
  }



  // TODO: Is this parallel? Can we make this batch-parallel or is it thread-safe?
  uintE insert(sequence<uintE> batch) {


    
    batch = pbbslib::sample_sort(batch, [&] (const uintE u, const uintE v) {
      return G->get_vertex(u).getOutDegree() > G->get_vertex(v).getOutDegree();
    });

    //Stores degrees and sorts in descending order
    sequence<uintE> deg = pbbs::sequence<uintE>(batch.size());
    par_for(0, batch.size(), [&] (size_t i) {
      deg[i] = deg[i] = G->get_vertex(batch[i]).getOutDegree();
    });

    //--------------------------Compute H-index--------------------------//

    pbbs::sequence<uintE> lowDegSum = pbbs::sequence<uintE>(0);
    pbbs::sequence<std::tuple<uintE, uintE>> highDegs;
    
   
    lowDegSum = lowDegD.slice(hindex + 1, hindex + batch.size() + 1);    



    //---------------------------ADDS TO C-------------------------------------//

    //Offsets prevents overriding current vertex entries

    //Resizes lowDegC and fills offsets
    // TODO: Why are you resizing each block in lowDegC? You should only resize if you
    // need the space. Figure out what space you need first, then do all of 
    // the necessary resizing.
    sequence<uintE> lowDegCOffset = pbbs::sequence<uintE>(lowDegC.size());
    par_for(0, lowDegC.size(), [&] (size_t i) {
      lowDegCOffset[i] = lowDegC[i]->size();
      lowDegC[i]->resize(batch.size() + lowDegC[i]->size(), UINT_E_MAX); //Fills with empty value which will be filtered later
    });

    // TODO: Why do you need to take out all of the elements in highDegC?
    //Resizes highDegC and fills offsets
    auto entries = highDegC.entries();
    //Will sparse_table slow this down?
    // TODO: Probably yes; if we can avoid it it'd be better.
    auto highDegCOffsets = make_sparse_table<uintE, uintE, hash_uintE>
      (entries.size(), std::make_tuple(UINT_E_MAX, UINT_E_MAX), hash_uintE());
    //Clears C, using entries, adds them back to C in parallel
    highDegC.clear();


    //Computes Prefix Sum (note scan_add_inplace is exclusive)
    //lowDeg
    pbbslib::scan_add_inplace(lowDegSum);
    par_for(0, batch.size() + 1, [&] (size_t i) {
      lowDegSum[i] += lowDegD[hindex + i + 1] + B.size(); //Compute inclusive and add size of B
    });

    

    //Creates the binary array described as M
    pbbs::sequence<bool> M = pbbs::sequence<bool>(lowDegSum.size() + 1);
  
    M[0] = (deg[B.size()] > hindex);
    par_for(0, lowDegSum.size(), [&] (size_t i) {
      if (lowDegSum[i] < deg.size()) {
        M[i + 1] = (deg[lowDegSum[i]] > hindex + 1 + i);
      }
      else {
        M[i + 1] = false;
      }
      
    });

          

    //Find index of first 0
    pbbs::sequence<size_t> indexM = pbbs::sequence<size_t>(M.size());
    par_for(0, M.size(), [&] (size_t i) {
      indexM[i] = !M[i] ? i : UINT_E_MAX;
    });

    

    
    
    //--------------------------ADD TO C--------------------------//
    //Subraction
    sequence<bool> difference = pbbs::sequence<bool>(batch.size() - 1);
    
    par_for(0, deg.size() - 1, [&] (size_t i) {
      difference[i] = (deg[i] != deg[i + 1]);  
    });

    sequence<size_t> indices = pbbs::sequence<size_t>(difference.size() + 1);

    par_for(0, difference.size() + 1, [&] (size_t i) {
      indices[i] = i;
    });

    auto f = [&] (size_t i) { return difference[i] || i == batch.size() - 1; };
    indices = pbbs::filter(indices, f);
    //Uses indices sequence to know the index of last entry for each clustered deg
    par_for(0, indices.size(), [&] (size_t i) {
      size_t start = (i == 0 ? 0 : indices[i - 1] + 1);
      size_t end = indices[i] + 1;

      pbbs::sequence<uintE> extra = batch.slice(start, end);
      
      if (deg[i] > threshold) {
        if (highDegC.find(deg[indices[i]], nullptr) == nullptr) {
          highDegC.insert(std::make_tuple(deg[indices[i]], &extra));
        }
        else {
          highDegC.find(deg[indices[i]], nullptr)->append(extra);
        }

        if (highDegD.find(deg[indices[i]], UINT_E_MAX) == UINT_E_MAX) {
          highDegD.insert(std::make_tuple(deg[indices[i]], extra.size()));
        }
        else {
          highDegD.change(deg[indices[i]], highDegD.find(deg[indices[i]], UINT_E_MAX) + extra.size()); 
        }
      }

      else {
        lowDegC[deg[indices[i]]]->append(extra);
        lowDegD[deg[indices[i]]] += extra.size();
      }
    });

    //--------------------------ADD TO C--------------------------//

    uintE hindexIncrease = pbbslib::reduce_min(indexM);
    if (hindexIncrease != UINT_E_MAX) {
      hindex += hindexIncrease;
    }
    else {
      hindex += batch.size();
    }

    return hindex;
  }

   
  uintE erase(uintE v) {
    //Implement later
    return hindex;
  }
	
  uintE change(uintE v) {
/*
    erase(v);
    insert(v);

    //Why only add elements to P in an update?
    //What is the point of maintaining the dictionary parallel to C? Can't we just add during the update?
    if (G->get_vertex(v).getOutDegree() >= 2 * hindex) {
      P.insert(std::make_tuple(v, pbbs::empty()));
    }
*/
    return hindex;
  }

  void addToC(uintE v) {
/*
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
        lowDegC[deg] = new pbbs::sequence<uintE>(G->n);
      }
      lowDegC[deg]->add(v);
    }
*/
  }

  void removeFromC(uintE v) { // Assumes v is in C
    /*
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
*/
  }

  bool containedInC(uintE deg) {
/*
    if (deg > threshold) {

      return highDegC.contains(deg);
    }
    else {
      if (lowDegC[deg] == nullptr) {
        return false;
      }
      return lowDegC[deg]->size != 0;
    }
   */ 
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
