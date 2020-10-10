

#include "hindex_dyn_arr.h"
#include "hindex_threshold.h"
#include "hindex_sparse_table.h"
#include "TriangleCounting.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "utils/generators/barabasi_albert.h"
#include "pbbslib/get_time.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <assert.h>

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
  std::cout << "### Params: -type = " << type << std::endl;
  std::cout << "### ------------------------------------" << endl;

  assert(P.getOption("-s"));
  assert(type < 2); //Valid option (will increase as there are more options)

  timer insertion;
  timer insertionTotal;
  timer deletion;
  timer deletionTotal;

  timer triangleTime;
  timer totalTime;
  
  
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
/*
  //Add all edges from static graph
  pbbs::sequence<uintE> degrees = pbbs::sequence<uintE>(GA.n);
  par_for(0, GA.n, [&] (size_t i) {
    degrees[i] = GA.get_vertex(i).degree;
  });

  uintE maxDeg = pbbslib::reduce_max(degrees);

  for (long idx = 0; idx < maxDeg; idx++) {

    pbbs::sequence<std::pair<uintE, uintE>> batch = pbbs::sequence<std::pair<uintE, uintE>>(GA.n);
    par_for(0, GA.n, [&] (size_t i) {
      if (GA.get_vertex(i).degree > 0) {
        batch[i] = std::make_pair(i, GA.get_vertex(i).getOutNeighbor(idx % GA.get_vertex(i).degree));
      }
      else {
        batch[i] = std::make_pair(UINT_E_MAX, UINT_E_MAX);
      }
    });

    triangle.addEdges(getEdges(batch));
    //if (idx % (maxDeg / 10) == 0) cout << "Triangles: " << triangle.total << endl;
  }

  cout << "Static Graph Complete" << endl;
*/
  totalTime.start();
  insertionTotal.start();

  //Add random edges
  for (int i = 0; i < 10; i++) {
    cout << "Batch " << (i + 1) << endl;
    //Random number of vertices between 10^2 to 10^3, each with 100 edges
    auto batch = barabasi_albert::generate_updates(100 + rand() % (1000 - 100), 100);

    triangleTime.start();
    insertion.start();
    triangle.addEdges(getEdges(batch));
    insertion.stop();
    triangleTime.stop();
  }
  insertionTotal.stop();
  cout << "Triangles: " << triangle.total << endl;

  deletionTotal.start();
  
  //Delete random edges
  for (int i = 0; i < 10; i++) {
    cout << "Batch " << (i + 11) << endl;
    //Random number of vertices between 10^2 to 10^3, each with 100 edges
    auto batch = barabasi_albert::generate_updates(100 + rand() % (10000 - 100), 100);

    triangleTime.start();
    deletion.start();
    triangle.removeEdges(getEdges(batch));
    deletion.stop();
    triangleTime.stop();
  }
  deletionTotal.stop();
  totalTime.stop();

  cout << "Triangles: " << triangle.total << endl;

  cout << "Actual Insertion Time: " << insertion.get_total() << endl;
  cout << "Total Insertion Time: " << insertionTotal.get_total() << endl;
  cout << "Actual Deletion Time: " << insertion.get_total() << endl;
  cout << "Total Deletion Time: " << insertionTotal.get_total() << endl;
  cout << "Actual Counting Time: " << triangleTime.get_total() << endl;
  cout << "Total Time: " << totalTime.get_total() << endl;

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
