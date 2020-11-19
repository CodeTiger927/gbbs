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

pbbs::sequence<std::pair<uintE, uintE>> getEdges(pbbs::sequence<std::pair<uintE, uintE>> edges) {

  auto edgesOrdered = pbbs::sequence<std::pair<uintE, uintE>>(edges.size(), [&] (size_t i) {
    if (edges[i].first > edges[i].second) return std::make_pair(edges[i].second, edges[i].first);
    else return edges[i];
  });

  //Get unique edges
  sequence<std::pair<uintE,uintE>> sortedE = merge_sort(edgesOrdered, [&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {
    if (a.first != b.first) return a.first < b.first;
    else return a.second < b.second;
  });

  auto eTemp = sequence<std::pair<uintE, uintE>>(sortedE.size(), [&] (size_t i) {
    if (i == 0) {
      return sortedE[i];
    }
    else {
      if (sortedE[i].first != sortedE[i - 1].first || sortedE[i].second != sortedE[i - 1].second) {
        return sortedE[i];
      }
      else return std::make_pair(UINT_E_MAX, UINT_E_MAX);
    }
  });    

  auto uniqueEdges = filter(eTemp, [&] (std::pair<uintE, uintE> i) { return (i.first != UINT_E_MAX && i.second != UINT_E_MAX) && (i.first != i.second); } );
 
  return uniqueEdges;
}

// -type 0 for HSetDyn_arr
// -type 1 for HSetThreshold

//e.g.: ./SubgraphCounting -s -rounds 1 -type 0 -size 10 "inputs/graph_test_3.txt"

template <class Graph>
double AppSubgraphCounting_runner(Graph& GA, commandLine P) {

  long type = static_cast<uintE>(P.getOptionLongValue("-type", 0));
  long size = static_cast<uintE>(P.getOptionLongValue("-size", 10));

  std::cout << "### Application: Subgraph Counting" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### Params: -type = " << type << " -size = " << size << std::endl;
  std::cout << "### ------------------------------------" << endl;

  assert(P.getOption("-s"));
  assert(type < 2); //Valid option (will increase as there are more options)

  // Parameters for barabsi_albert
  // Parameter 1  
  uintE barabasi_albert_parameter1 = 10;
  // Parameter 2
  uintE barabasi_albert_parameter2 = 5;

  timer insertion;
  timer insertionTotal;
  timer deletion;
  timer deletionTotal;

  timer staticTime;
  timer triangleTime;
  timer totalTime;
  
  
  auto _dynG = createEmptyDynamicSymmetricGraph<dynamic_symmetric_vertex, pbbs::empty>();
  HSet* h;

  if (type == 0) {
    h = new HSetDynArr(&_dynG);
    std::cout << "DYN_ARR VERSION\n" << std::endl;
  }
  //else if (type == 1) {
    //h = new HSetThreshold(&_dynG, GA.n);

    //std::cout << "THRESHOLD VERSION\n" << std::endl;
  //}

  TriangleCounting triangle = TriangleCounting(h);
  // Temporary for testing purposes. Initialize this to max node.
  triangle.initialize(GA.n + 1000);

  
  sequence<uintE> sizes = sequence<uintE>(GA.n,[&](size_t i) {return GA.get_vertex(i).degree;});
  pbbslib::scan_add_inplace(sizes);
  pbbslib::dyn_arr<std::pair<uintE,uintE>> insertBatchtmp = pbbslib::dyn_arr<std::pair<uintE,uintE>>(GA.m);
  par_for(0,GA.n,[&](size_t i) {
    par_for(0,GA.get_vertex(i).degree,[&](size_t j) {
      uintE n = GA.get_vertex(i).getOutNeighbor(j);
      if(i > n) {
        insertBatchtmp.A[sizes[i] + j] = std::make_pair(UINT_E_MAX,UINT_E_MAX);
      }else{
        insertBatchtmp.A[sizes[i] + j] = std::make_pair(i,n);
      }
    });
  });
  insertBatchtmp.size = insertBatchtmp.capacity;
  sequence<std::pair<uintE,uintE>> insertBatch = filter(insertBatchtmp.to_seq(),[&](std::pair<uintE,uintE> p) {return (p.first != UINT_E_MAX || p.second != UINT_E_MAX);});

  totalTime.start();

  // Add all edges from static graph
  staticTime.start();
  triangle.addEdges(insertBatch);
  staticTime.stop();
  

  std::cout << "Initial Triangle Count: " << triangle.total << std::endl;

  insertionTotal.start();

  //Add random edges
  for (int i = 0; i < size; i++) {
    cout << "Batch " << (i + 1) << endl;
    //Random number of vertices between 10^2 to 10^3, each with 100 edges
    auto batch = getEdges(barabasi_albert::generate_updates(barabasi_albert_parameter1, barabasi_albert_parameter2));

    triangleTime.start();
    insertion.start();
    triangle.addEdges(batch);
    insertion.stop();
    triangleTime.stop();
    batch.clear();
    // cout << "Triangles: " << triangle.total << endl;
  }
  insertionTotal.stop();
  cout << "Triangles: " << triangle.total << endl;

  deletionTotal.start();
  
  //Delete random edges
  for (int i = 0; i < size;i++) {
    cout << "Batch " << (i + size + 1) << endl;
    //Random number of vertices between 10^2 to 10^3, each with 100 edges
    auto batch = getEdges(barabasi_albert::generate_updates(barabasi_albert_parameter1, barabasi_albert_parameter2));

    triangleTime.start();
    deletion.start();
    triangle.removeEdges(batch);
    deletion.stop();
    triangleTime.stop();
    batch.clear();
  }
  deletionTotal.stop();
  totalTime.stop();

  cout << "Triangles: " << triangle.total << endl;
  cout << "Inserting Static Graph time: " << staticTime.get_total() << endl;
  cout << "Actual Insertion Time: " << insertion.get_total() << endl;
  cout << "Total Insertion Time: " << insertionTotal.get_total() << endl;
  cout << "Actual Deletion Time: " << insertion.get_total() << endl;
  cout << "Total Deletion Time: " << insertionTotal.get_total() << endl;
  cout << "Actual Counting Time: " << triangleTime.get_total() << endl;
  cout << "Total Time: " << totalTime.get_total() << endl;

  h->del();

  par_for(0, h->G->n, [&] (size_t i) {
    if (h->G->existVertices.A[i]) pbbslib::free_array(h->G->v_data.A[i].neighbors.table);
  });

  pbbslib::free_array(h->G->v_data.A);
  pbbslib::free_array(h->G->existVertices.A);

  triangle.del();

  return 0;
}

/* RANDOM EDGE RUNNER
template <class Graph>
double AppSubgraphCounting_runner(Graph& GA, commandLine P) {
  long type = static_cast<uintE>(P.getOptionLongValue("-type", 0));
  long size = static_cast<uintE>(P.getOptionLongValue("-size", 10));
  std::cout << "### Application: Subgraph Counting" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### Params: -type = " << type << std::endl;
  std::cout << "### ------------------------------------" << endl;
  assert(P.getOption("-s"));
  assert(type < 3); //Valid option (will increase as there are more options)
  timer clock;
  //auto _dynG = dynamifyDSG<dynamic_symmetric_vertex, pbbs::empty, Graph>(GA); //Dynamify inside main
  auto _dynG = createEmptyDynamicSymmetricGraph<dynamic_symmetric_vertex, pbbs::empty>();
  HSet* h;
  if (type == 0) {
    h = new HSetDynArr(&_dynG);
    std::cout << "DYN_ARR VERSION\n" << std::endl;
  }
  else if (type == 1) {
    h = new HSetThreshold(&_dynG, GA.n);
    std::cout << "THRESHOLD VERSION\n" << std::endl;
  }
  TriangleCounting triangle = TriangleCounting(h);
  clock.start();
  pbbs::sequence<pbbs::sequence<std::pair<uintE, uintE>>> batches = pbbs::sequence<pbbs::sequence<std::pair<uintE, uintE>>>(size);
  for (int i = 0; i < size; i++) {
    std::cout << "\n-----CHANGE " << (i + 1) << "-----" << std::endl;
    batches[i] = barabasi_albert::generate_updates(rand() % 1000 + 50, rand() % 100 + 5);
    //batches[i] = barabasi_albert::generate_updates(rand() % 10, rand() % 5); //Smaller dataset easier to debug
    //h->insertEdges(getEdges(batches[i]));
    triangle.addEdges(getEdges(batches[i]));
    int t = 0;
    for (int u = 0; u < h->G->n; u++) {
      for (int v = u + 1; v < h->G->n; v++) {
        for (int w = v + 1; w < h->G->n; w++) {
          if (h->G->existEdge(u, v) && h->G->existEdge(u, w) && h->G->existEdge(v, w)) t++;
        }
      }
    }
    std::cout << "Trianges: " << triangle.total << std::endl;
    std::cout << "h-index: " << h->hindex << std::endl;
    assert(t == triangle.total);
  }
  
  for (int i = size - 1; i >= 0; i--) {
    std::cout << "\n-----CHANGE " << (2 * size - i) << "-----" << std::endl;
    //h->eraseEdges(getEdges(batches[i]));
    triangle.removeEdges(getEdges(batches[i]));
    int t = 0;
    for (int u = 0; u < h->G->n; u++) {
      for (int v = u + 1; v < h->G->n; v++) {
        for (int w = v + 1; w < h->G->n; w++) {
          if (h->G->existEdge(u, v) && h->G->existEdge(u, w) && h->G->existEdge(v, w)) t++;
        }
      }
    }
    std::cout << "Trianges: " << triangle.total << std::endl;
    std::cout << "h-index: " << h->hindex << std::endl;
    assert(t == triangle.total);
  }
  
  std::cout << "\n\n------------------------------------" << endl;
  std::cout << "FINAL GRAPH SIZE: " << h->G->n << ", " << h->G->m << endl;
  std::cout << "RUN TIME: " << clock.stop() << std::endl;
  return 0;
}

*/
generate_symmetric_main(AppSubgraphCounting_runner, false);
