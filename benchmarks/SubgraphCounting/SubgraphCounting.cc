#include "hindex_dyn_arr.h"
//#include "hindex_threshold.h"
#include "TriangleCounting.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "utils/generators/barabasi_albert.h"
#include "pbbslib/get_time.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <assert.h>
#include <vector>


// Map the dense components to spread out throughout the graph
pbbs::sequence<std::pair<uintE, uintE>> sparcifyEdges
  (pbbs::sequence<std::pair<uintE, uintE>> edges,uintE l,uintE h,
    pbbs::random& rand) {
  auto sparseEs = make_sparse_table<uintE,uintE,hash_uintE>(2 * edges.size(),
    std::make_tuple(UINT_E_MAX,UINT_E_MAX),hash_uintE());
  par_for(0,edges.size(),[&](size_t i) {
    uintE u = edges[i].first;
    uintE v = edges[i].second;
    if(sparseEs.find(u,UINT_E_MAX) == UINT_E_MAX) {
      sparseEs.insert(std::make_tuple(u,rand[2 * i] % (h - l + 1) + l));
    }
    if(sparseEs.find(v,UINT_E_MAX) == UINT_E_MAX) {
      sparseEs.insert(std::make_tuple(v,rand[2 * i + 1] % (h - l + 1) + l));
    }
  });
  pbbs::sequence<std::pair<uintE, uintE>> res(edges.size(),
    [&](size_t i) {return std::make_pair(
      sparseEs.find(edges[i].first,0),sparseEs.find(edges[i].second,0));});
  sparseEs.del();
  rand = rand.next();
  return res;
}

pbbs::sequence<std::pair<uintE, uintE>> getEdges
  (pbbs::sequence<std::pair<uintE, uintE>> edges,uintE N) {
  pbbs::random rand = pbbs::random();
  edges = sparcifyEdges(edges,0,N - 1,rand);

  auto edgesOrdered = pbbs::sequence<std::pair<uintE, uintE>>
    (edges.size(), [&] (size_t i) {
      if (edges[i].first > edges[i].second) {
        return std::make_pair(edges[i].second, edges[i].first);
      }
      else return edges[i];
    }
  );

  //Get unique edges
  sequence<std::pair<uintE,uintE>> sortedE = merge_sort(edgesOrdered, 
    [&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {
      if (a.first != b.first) return a.first < b.first;
      else return a.second < b.second;
    }
  );

  auto eTemp = sequence<std::pair<uintE, uintE>>(sortedE.size(), 
    [&] (size_t i) {
      if (i == 0) {
        return sortedE[i];
      }
      else {
        if (sortedE[i].first != sortedE[i - 1].first || 
            sortedE[i].second != sortedE[i - 1].second) {

          return sortedE[i];
        }
        else return std::make_pair(UINT_E_MAX, UINT_E_MAX);
      }
    }
  );    

  auto uniqueEdges = filter(eTemp, [&] (std::pair<uintE, uintE> i) { 
    return (i.first != UINT_E_MAX && i.second != UINT_E_MAX) && (i.first != i.second);
  });
 
  return uniqueEdges;
}


// -type 0 for HSetDyn_arr
// -type 1 for HSetThreshold

// ./SubgraphCounting -s -rounds 1 -type 0 -size 10 -nodes 100 -edges 10
//   "inputs/graph_test_3.txt"

template <class Graph>
double AppSubgraphCounting_runner(Graph& GA, commandLine P) {

  long type = static_cast<uintE>(P.getOptionLongValue("-type", 0));
  long size = static_cast<uintE>(P.getOptionLongValue("-size", 10));
  //Barabasi_albert parameter - nodes
  uintE nodes = static_cast<uintE>(P.getOptionLongValue("-nodes", GA.n));
  //Barabasi_albert parameter - edges
  uintE edges = static_cast<uintE>(P.getOptionLongValue("-edges", GA.m / GA.n));
  // Mode activated for testing to more easily store output
  bool scriptMode = static_cast<bool>(P.getOptionLongValue("-scriptMode",false));


  if(!scriptMode) {
    std::cout << "### Application: Subgraph Counting" << std::endl;
    std::cout << "### Graph: " << P.getArgument(0) << std::endl;
    std::cout << "### Threads: " << num_workers() << std::endl;
    std::cout << "### n: " << GA.n << std::endl;
    std::cout << "### m: " << GA.m << std::endl;
    std::cout << "### Params: -type = " << type << " -size = " << size << std::endl;
    std::cout << "### ------------------------------------" << endl;
  }
  assert(P.getOption("-s"));
  assert(type < 2); //Valid option (will increase as there are more options)
  
  timer insertion;
  timer insertionTotal;
  timer deletion;
  timer deletionTotal;

  timer staticTime;
  timer triangleTime;
  timer totalTime;
  
  auto _dynG = createEmptyDynamicSymmetricGraph
    <dynamic_symmetric_vertex, pbbs::empty>();

  HSet* h;

  if (type == 0) {
    h = new HSetDynArr(&_dynG);
    std::cout << "DYN_ARR VERSION\n" << std::endl;
  }
  //else if (type == 1) { //Temporarily disable threshold (needs work)
    //h = new HSetThreshold(&_dynG, GA.n);

    //std::cout << "THRESHOLD VERSION\n" << std::endl;
  //}

  TriangleCounting triangle = TriangleCounting(h, false);
  // Temporary for testing purposes. Initialize this to max node.
  triangle.initialize(GA.n + 1000);

  
  sequence<uintE> sizes = sequence<uintE>(GA.n,[&](size_t i) {
    return GA.get_vertex(i).degree;
  });

  pbbslib::scan_add_inplace(sizes);
  sequence<std::pair<uintE,uintE>> insertBatchtmp = 
    sequence<std::pair<uintE,uintE>>(GA.m);


  par_for(0,GA.n,[&](size_t i) {
    par_for(0,GA.get_vertex(i).degree,[&](size_t j) {

      uintE n = GA.get_vertex(i).getOutNeighbor(j);
      if(i > n) {
        insertBatchtmp[sizes[i] + j] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
      }else{
        insertBatchtmp[sizes[i] + j] = std::make_pair(i,n);
      }
    });
  });


  sequence<std::pair<uintE,uintE>> insertBatch = 
    filter(insertBatchtmp, [&](std::pair<uintE,uintE> p) {
      return (p.first != UINT_E_MAX || p.second != UINT_E_MAX);
    }
  );

  totalTime.start();

  // Add all edges from static graph
  staticTime.start();
  triangle.addEdges(insertBatch);
  staticTime.stop();
  uintE staticHIndex = triangle.getHIndex();
  uintE staticTriangle = triangle.total;
  if(!scriptMode) std::cout << "Initial Triangle Count: " << triangle.total << std::endl;
  insertionTotal.start();

  //Add random edges
  for (int i = 0; i < size; i++) {
    if (i % 10 == 0 && !scriptMode) cout << "Batch " << (i + 1) << endl;
    
    auto batch = getEdges(
      barabasi_albert::generate_updates(nodes, edges),GA.n
    );

    triangleTime.start();
    insertion.start();
    triangle.addEdges(batch);
    insertion.stop();
    triangleTime.stop();
    batch.clear();
    // cout << "Triangles: " << triangle.total << endl;
  }
  insertionTotal.stop();
  uintE insertionHIndex = triangle.getHIndex();
  uintE insertionTriangle = triangle.total;

  if(!scriptMode) cout << "Triangles: " << triangle.total << endl;

  deletionTotal.start();
  
  //Delete random edges
  for (int i = 0; i < size;i++) {
    if (i % 10 == 0 && !scriptMode) cout << "Batch " << (i + size + 1) << endl;
    //Random number of vertices between 10^2 to 10^3, each with 100 edges
    auto batch = getEdges(
      barabasi_albert::generate_updates(nodes, edges),GA.n
    );

    triangleTime.start();
    deletion.start();
    triangle.removeEdges(batch);
    deletion.stop();
    triangleTime.stop();
    batch.clear();
  }
  deletionTotal.stop();
  uintE deletionHIndex = triangle.getHIndex();
  uintE deletionTriangle = triangle.total;
  totalTime.stop();

  if(!scriptMode) {
    cout << "Triangles: " << triangle.total << endl;
    cout << "Inserting Static Graph time: " << staticTime.get_total() << endl;
    cout << "Actual Insertion Time: " << insertion.get_total() << endl;
    cout << "Total Insertion Time: " << insertionTotal.get_total() << endl;
    cout << "Actual Deletion Time: " << insertion.get_total() << endl;
    cout << "Total Deletion Time: " << insertionTotal.get_total() << endl;
    cout << "Actual Counting Time: " << triangleTime.get_total() << endl;
    cout << "Total Time: " << totalTime.get_total() << endl;
  }

  h->del();

  par_for(0, h->G->n, [&] (size_t i) {
    if (h->G->existVertices.A[i]) {
      pbbslib::free_array(h->G->v_data.A[i].neighbors.table);
    }
  });

  pbbslib::free_array(h->G->v_data.A);
  pbbslib::free_array(h->G->existVertices.A);

  triangle.del();
  
  if(scriptMode) {
      cout << type << ", " << size << ", " << nodes << ", " << edges << ", " 
    << staticTime.get_total() << ", " << insertionTotal.get_total() << ", "
    << deletionTotal.get_total() << ", " << totalTime.get_total() << ", "
    << staticTriangle << ", " << insertionTriangle << ", " << deletionTriangle
    << endl;
  }

  return 0;
}

generate_symmetric_main(AppSubgraphCounting_runner, false);
