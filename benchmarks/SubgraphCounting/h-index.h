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

// For hashing sparse table

//Already defined
//struct hash_uintE {
  //inline size_t operator () (const uintE & a) {return pbbs::hash64_2(a);}
//};

template <class Graph>
struct HSet {

  Graph* G; //Graph

  sparse_table<uintE, pbbs::empty, hash_uintE> H; //Set of elements such that deg(x) >= |H|
  size_t hindex; //|H|

  sparse_table<uintE, pbbs::empty, hash_uintE> P;

  pbbs::sequence<uintE> B; //Set of elements such that deg(x) == |H|

  pbbs::sequence<pbbs::sequence<uintE>*> lowDegC; //Map of elements not in H fed by degree up to a threshold
  sparse_table<uintE, pbbs::sequence<uintE>*, hash_uintE> highDegC; // Map of elements not in H sorted by degree greater than the threshold
  uintE threshold;

  pbbs::sequence<uintE> lowDegD;
  sparse_table<uintE, uintE, hash_uintE> highDegD;

  //Doesn't have to be symmetric vertex
  //Specify templatized graph
  //a - threshold for low and high degree vertices
  HSet(Graph& _G) {

    G = &_G;
    threshold = G->n;

    
    H = make_sparse_table<uintE, pbbs::empty, hash_uintE> //|H| has to be between m/n and sqrt(2m)
      (std::max(G->m/G->n, (size_t) sqrt(2 * G->m)), std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());
    hindex = 0;

    P = make_sparse_table<uintE, pbbs::empty, hash_uintE>
      (std::max(G->m/G->n, (size_t) sqrt(2 * G->m)), std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());

    B = pbbs::sequence<uintE>(0);

    lowDegC = sequence<pbbs::sequence<uintE>*>(threshold + 1); //possible degrees range from 0 ... a
    // TODO: You can reduce the search space of highDegC by taking G->n - threshold
    // It won't make a huge difference, but maybe it'll be better for the hash
    // function to distribute these values evenly, since we don't have any
    // values in ehre from 0 .. threshold
    highDegC = make_sparse_table<uintE, pbbs::sequence<uintE>*, hash_uintE>
      (G->n, std::make_tuple(UINT_E_MAX, nullptr), hash_uintE());

    par_for(0, lowDegC.size(), [&] (size_t i) { lowDegC[i] = new sequence<uintE>(0); });

    lowDegD = sequence<uintE>(threshold + 1);
    par_for(0, lowDegD.size(), [&] (size_t i) { lowDegD[i] = 0; });
    highDegD = make_sparse_table<uintE, uintE, hash_uintE>
      (G->n, std::make_tuple(UINT_E_MAX, UINT_E_MAX), hash_uintE());
  }

//--------------------------INSERT--------------------------//
  uintE insert(sequence<uintE> batch) {
    
    //TODO: Consider vertices whose old degree and new degree are both above hindex
    //TODO: Update sets B and H

    // TODO: Why do you need to do batch.slice()? Versus just batch -- and same with .rslice()
    // TODO: You shouldn't be assigning to batch like this -- you're not properly getting rid
    // of the memory in batch and you're doing some weird copy assignment here.
    // Instead, you should have an out array and properly delete batch, and use the sorted
    // out array.
    batch = integer_sort(batch.slice(), [&] (uintE v) { return G->get_vertex(v).getOutDegree(); }).rslice();
    //Stores degrees and sorts in descending order
    sequence<uintE> deg = pbbs::sequence<uintE>(batch.size());
    par_for(0, batch.size(), [&] (size_t i) {
      // TODO: What's with this double assignment? What are you trying to do?
      deg[i] = deg[i] = G->get_vertex(batch[i]).getOutDegree();
    });

    //--------------------------Compute H-index--------------------------//

    pbbs::sequence<uintE> lowDegSum = pbbs::sequence<uintE>(0);
    pbbs::sequence<std::tuple<uintE, uintE>> highDegs = pbbs::sequence<std::tuple<uintE, uintE>>(0);

    
    // TODO: What is low deg d and why do you need it? How is it different from C?
    // TODO: It might be better if you just had the simpler low degree code here, instead of
    // trying to do this fancy thing with lowDeg vs highDeg -- especially since we haven't
    // even tested this and ensured that it works
    if (hindex + batch.size() <= threshold) { //hindex to hindex + batch.size() completely in lowDeg
      lowDegSum = lowDegD.slice(hindex + 1, hindex + batch.size() + 1);
    }
    else {
      if (hindex <= threshold) { //Range in both
        lowDegSum = lowDegD.slice(hindex + 1, lowDegD.size());

        auto inRange = [&] (std::tuple<uintE, uintE> deg) { return std::get<0>(deg) > threshold && std::get<0>(deg) <= hindex + batch.size(); };
        highDegs = pbbs::filter(highDegD.entries(), inRange);
    
      }
      else { //hindex to hindex + batch.size() completely in lowDeg
        auto inRange = [&] (std::tuple<uintE, uintE> deg) { return std::get<0>(deg) > hindex && std::get<0>(deg) <= hindex + batch.size(); };
        highDegs = pbbs::filter(highDegD.entries(), inRange);
      }
    }
    /*
    highDegs = pbbslib::sample_sort(highDegs, [&] (const std::tuple<uintE, uintE> u, const std::tuple<uintE, uintE> v) {
      return std::get<0>(u) < std::get<0>(v);
    });
    */
    highDegs = integer_sort(highDegs.slice(), [&] (std::tuple<uintE, uintE> v) { return std::get<0>(v); });
    //Computes Prefix Sum (note scan_add_inplace is exclusive)
    //lowDeg (even if lowDeg is empty it still checks whether or not it can empty B out)
    pbbslib::scan_add_inplace(lowDegSum);
    par_for(0, batch.size() + 1, [&] (size_t i) {
      lowDegSum[i] += lowDegD[hindex + i + 1]; //Compute inclusive sunm
    });

    // TODO: What is M?
    //Creates the binary array described as M
    pbbs::sequence<bool> M = pbbs::sequence<bool>(lowDegSum.size() + 1);

    M[0] = (deg[B.size()] > hindex);
    par_for(0, lowDegSum.size(), [&] (size_t i) {
      // TODO: Why do you have the +1 here? Just make your parallel for in the correct range.
      // TODO: Don't need to do an if else like this; just use the full boolean expression --
      // lowDegSum[i] < deg.size() && whatever.
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
    updateC(batch, deg);

    //--------------------------RETURN HINDEX--------------------------//
    uintE hindexIncrease = pbbslib::reduce_min(indexM);
    // TODO: What you're doing here doesn't make sense to me. Explain it better. Also, where
    // do you update H? What about updating D? Why was C never used?


    if (hindex == UINT_E_MAX) { //Check highDegD because new hindex not in lowDeg

      pbbs::sequence<uintE> highDegSum;
      if (highDegs.size() == 0) {
        highDegs = pbbs::sequence<std::tuple<uintE, uintE>>(1);
        highDegs[0] = std::make_tuple(0,0);

        highDegSum = pbbs::sequence<uintE>(1);
        highDegSum[0] = 0;
      }
      else {
        highDegSum = pbbs::sequence<uintE>(highDegs.size());
        par_for(0, highDegs.size(), [&] (size_t i) {
          highDegSum[i] = std::get<1>(highDegs[i]);
        });

        //Inclusive Prefix sum
        pbbslib::scan_add_inplace(highDegSum);
        par_for(0, highDegs.size(), [&] (size_t i) {
          highDegSum[i] += std::get<1>(highDegs[i]) + lowDegSum[lowDegSum.size() - 1]; //Adds on last element of lowDegSum
        });
      }
      

      
      auto highM = pbbs::sequence<uintE>(highDegs.size());

      par_for(0, highM.size() - 1, [&] (size_t i) {
        //Interval:        deg1              to            deg 2
        //          std::get<0>(highDegs[i]) to std::get<0>(highDegs[i + 1])

        if (highDegSum[i] >= deg.size() || deg[highDegSum[i]] < std::get<0>(highDegs[i])) { //New hindex is below this interval
          highM[i] = std::get<0>(highDegs[i]);
        }
        else if (deg[highDegSum[i]] > std::get<0>(highDegs[i + 1])) { //New hindex is above this interval
          highM[i] = UINT_E_MAX; //Flag, equivalent to true in lowDeg M
        }
        else { //New hindex in this interval
          highM[i] = deg[highDegSum[i]];
        }
      });
                            
      //Last interval: final degree in highDegs to hindex + batch.size()
      if (highDegSum[highDegs.size() - 1] >= deg.size() || deg[highDegSum[highDegs.size() - 1]] < std::get<0>(highDegs[highDegs.size() - 1])) { //New hindex is below this interval
        highM[highDegs.size() - 1] = std::get<0>(highDegs[highDegs.size() - 1]);
      }
      //hindex cannot be above this range

      else { //New hindex in this interval
        highM[highDegs.size() - 1] = deg[highDegSum[highDegs.size() - 1]];
      }

      //Finding first non number that is not UINT_E_MAX
      auto indexHighM = pbbs::sequence<uintE>(highM.size());
      par_for(0, highM.size(), [&] (size_t i) {
        indexHighM[i] = highM[i] != UINT_E_MAX ? i : UINT_E_MAX;
      });

      auto index = pbbslib::reduce_min(indexHighM);

      if (index == UINT_E_MAX) {
        hindex += batch.size();
      }
      else {
        hindex = highM[index];
      }

    }

    else {
      hindex += hindexIncrease;
    }

    
    return hindex;
  }

//--------------------------ERASE--------------------------//   
  uintE erase(sequence<uintE> batch) {
    return hindex;
  }

  //--------------------------ADD TO C--------------------------//
  void updateC(pbbs::sequence<uintE> batch, pbbs::sequence<uintE> deg) {
    
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
        // TODO: Why do you need C and D to be separate? What's the point? Also, where do you use C?
        lowDegC[deg[indices[i]]]->append(extra);
        lowDegD[deg[indices[i]]] += extra.size();
      }
    });
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
