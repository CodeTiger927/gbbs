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

  //pbbs::sequence<uintE> B; //Set of elements such that deg(x) == |H| Don't think I need this for the batch insert

  pbbs::sequence<pbbs::sequence<uintE>*> lowDegC; //Map of elements not in H sorted by degree up to a threshold
  sparse_table<uintE, pbbs::sequence<uintE>*, hash_uintE> highDegC; // Map of elements not in H sorted by degree greater than the threshold
  uintE threshold;


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

    //B = pbbs::sequence<uintE>(G->n);

    // C = make_sparse_table<uintE, sparse_table<uintE, pbbs::empty, hash_uintE>, hash_uintE>
      // (n, std::make_tuple(UINT_E_MAX, empty), hash_uintE());

    lowDegC = sequence<pbbs::sequence<uintE>*>(a); //possible degrees range from 0 ... a - 1
    // TODO: You can reduce the search space of highDegC by taking G->n - threshold
    // It won't make a huge difference, but maybe it'll be better for the hash
    // function to distribute these values evenly, since we don't have any
    // values in ehre from 0 .. threshold
    highDegC = make_sparse_table<uintE, pbbs::sequence<uintE>*, hash_uintE>
      (G->n, std::make_tuple(UINT_E_MAX, nullptr), hash_uintE());

    par_for(0, lowDegC.size(), [&] (size_t i) { lowDegC[i] = new sequence<uintE>(0); });

  }

  // TODO: Is this parallel? Can we make this batch-parallel or is it thread-safe?
  uintE insert(sequence<uintE> batch) {

    //Adding stuff to C
    //This might be really inefficient because it allocates space for each vertex in every entry for C
    //Then it filters out the ones not used (which are most of them)


    sequence<uintE> deg = pbbs::sequence<uintE>(batch.size());


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

    //Resizes arrays and then adds them back to C
    par_for(0, entries.size(), [&] (size_t i) {
      highDegCOffsets.insert(std::make_tuple(std::get<0>(entries[i]), std::get<1>(entries[i])->size()));
      std::get<1>(entries[i])->resize(batch.size() + std::get<1>(entries[i])->size(), UINT_E_MAX); //Fills with empty values which will be filtered later
      highDegC.insert(std::make_tuple(std::get<0>(entries[i]), std::get<1>(entries[i])));
    });


    //Adds vertices to C
    par_for(0, batch.size(), [&] (size_t i) {
    //for (size_t i = 0; i < batch.size(); i++) {

      deg[i] = G->get_vertex(batch[i]).getOutDegree();
      if (deg[i] < threshold) {
        (*lowDegC[deg[i]])[i + lowDegCOffset[deg[i]]] = batch[i];
      }
      else {

        //Degree entry in C
        if (highDegC.find(deg[i], nullptr) != nullptr) {

          //Degree already had entry in C
          if (highDegCOffsets.find(deg[i], UINT_E_MAX) != UINT_E_MAX) {
            (*highDegC.find(deg[i], nullptr))[i + highDegCOffsets.find(deg[i], UINT_E_MAX)] = batch[i];
          }
          //Degree entry was just added in this batch
          else  {
            (*highDegC.find(deg[i], nullptr))[i] = batch[i];
          }
        }
        
        //Degree entry not in C, needs to add
        else {

          //I don't know if this will cause problems in parallel
          sequence<uintE>* s = new pbbs::sequence<uintE>(0);
          s->resize(batch.size(), UINT_E_MAX);
          highDegC.insert(std::make_tuple(deg[i], s));
          //highDegC.find(deg[i], nullptr)->resize(batch.size(), UINT_E_MAX);
          (*highDegC.find(deg[i], nullptr))[i] = batch[i];
        }
      }
    });
    //}
    

    //Filters out empty spaces (removes gaps)
    //Offsets re-used because we need the new sizes
    auto f = [&] (uintE i) { return i != UINT_E_MAX; };

    par_for(0, lowDegC.size(), [&] (size_t idx) {
      auto filtered = pbbslib::filter(*lowDegC[idx], f);
      lowDegC[idx] = &filtered;
      lowDegCOffset[idx] = lowDegC[idx]->size();
    });

    entries = highDegC.entries();
    entries = pbbslib::sample_sort(entries, [&](const std::tuple<uintE, pbbs::sequence<uintE>*> deg1, const std::tuple<uintE, pbbs::sequence<uintE>*> deg2) {
      return std::get<0>(deg1) > std::get<0>(deg2); //Greater than gives degrees from greatest to least
    }, false); //Don't think I need it to be stable

    //I tried "highDegCOffsets.clear()" but for some reason it would take a really long time to execute line 164 (calling insert)
    sequence<std::tuple<uintE, uintE>> highDegCSizes = sequence<std::tuple<uintE, uintE>>(entries.size());
    
    highDegC.clear(); //Clears C so they can be readded
    par_for(0, entries.size(), [&] (size_t idx) {
    //for (int idx = 0; idx < entries.size(); idx++) {
      auto filtered = pbbs::filter(*std::get<1>(entries[idx]), f);
      highDegC.insert(std::make_tuple(std::get<0>(entries[idx]), &filtered));
      highDegCSizes[idx] = std::make_tuple(std::get<0>(entries[idx]), filtered.size());
      
    //}
    });

    
    //-------------------------Computes H-index--------------------------------//
    
    //auto highDegCSizes = highDegCOffsets.entries();
    
    auto fx = [&] (std::tuple<uintE, uintE> T) { return std::get<0>(T) > hindex; };
    highDegCSizes = pbbs::filter(highDegCSizes, fx); //Filters out all degrees that are less than H
    
    
    
    
    sequence<uintE> prefix_sum = sequence<uintE>(highDegCSizes.size());

    par_for(0, prefix_sum.size(), [&] (size_t i){
      prefix_sum[i] = std::get<1>(highDegCSizes[i]);
    });
    //Prefix sum so each degree knows how many vertices have degree greater than it
    
    auto temp = sequence<uintE>(prefix_sum);

    cout << "SIZE   " << highDegC.find(3, nullptr)->size() << endl; //C is fine right here
    pbbslib::scan_add_inplace(prefix_sum); //Computes exclusive prefix sum
    cout << "SIZE   " << highDegC.find(3, nullptr)->size() << endl; //But it gets scrambled heregf
    
    par_for(0, prefix_sum.size(), [&] (size_t i) {
      //Subtracts degree so all nonnegative numbers are potential h-indices
      prefix_sum[i] += std::get<1>(highDegCSizes[i]) - std::get<0>(highDegCSizes[i]); //Adds the size of it's own entry because prefix sum was exclusive
    });
    
    auto hFilter = [&] (uintE size) { return size >= 0; };
    auto hindices = pbbs::filter(prefix_sum, hFilter); //Possible h-indices
    
    
    //H-index is less than threshold (so search in lowDegC)
    if (hindices.size() == 0) {

      prefix_sum = lowDegCOffset.rslice(lowDegCOffset.size(), hindex);
      pbbslib::scan_add_inplace(prefix_sum);
      hindices = pbbs::filter(prefix_sum, f);

      par_for(0, prefix_sum.size(), [&] (size_t i) {
        //Subtracts degree so all nonnegative numbers are potential h-indices
        prefix_sum[i] += lowDegCOffset[i] - i; //Adds the size of it's own entry because prefix sum was exclusive
      });
      hindices = pbbs::filter(prefix_sum, hFilter);

      hindex = pbbslib::binary_search(prefix_sum, hindices[0], std::less<size_t>());
      
      H = make_sparse_table<uintE, pbbs::empty, hash_uintE>
        (2 * hindex, std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());

      //THIS DOESN'T WORK YET (supposed to update H)
      par_for(0, lowDegC.size(), [&] (size_t i) {
        if (i >= hindex) {

          par_for(0, lowDegC[i]->size(), [&] (size_t v) {
            H.insert(std::make_tuple((*lowDegC[i])[v], pbbs::empty()));
          });
        }
      });


    }
    else {
      size_t idx = pbbslib::binary_search(prefix_sum, hindices[0], std::less<size_t>());
      hindex = std::get<0>(highDegCSizes[idx]);
      H = make_sparse_table<uintE, pbbs::empty, hash_uintE>
        (2 * hindex, std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());
    }
    //THIS DOESN'T WORK YET (supposed to update H)
    par_for(0, entries.size(), [&] (size_t i) {
    //for(size_t i = 0; i < entries.size(); i++) {
      if (std::get<0>(entries[i]) >= hindex) {
        par_for(0, std::get<1>(entries[i])->size(), [&] (size_t v) {
        //for (int v = 0; v < std::get<1>(entries[i])->size(); v++) {

          H.insert(std::make_tuple((*std::get<1>(entries[i]))[v], pbbs::empty())); //SOMETHING HERE IS WRONG

        });
        //}
      }
    //}
    });
    
   
    return hindex;

  }


   
  uintE erase(uintE v) {
/*
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
*/
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
