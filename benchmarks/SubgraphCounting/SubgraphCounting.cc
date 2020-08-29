#include "hindex_dyn_arr.h"
#include "hindex_threshold.h"
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

  auto uniqueEdges = filter(eTemp, [&] (std::pair<uintE, uintE> i) { return (i.first != UINT_E_MAX && i.second != UINT_E_MAX); } );

  return uniqueEdges;
}

// -type 0 for HSetDyn_arr
// -type 1 for HSetThreshold

//e.g.: ./SubgraphCounting -s -rounds 1 -type 0 -size 10 "inputs/graph_test_3.txt"

template <class Graph>
double AppSubgraphCounting_runner(Graph& GA, commandLine P) {
  uintE type = static_cast<uintE>(P.getOptionLongValue("-type", 0));
  uintE size = static_cast<uintE>(P.getOptionLongValue("-size", 10));

  std::cout << "### Application: Subgraph Counting" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### Params: -type = " << type << std::endl;
  std::cout << "### ------------------------------------" << endl;

  assert(P.getOption("-s"));
  assert(type < 2); //Valid option (will increase as there are more options)

  timer clock;

  auto _dynG = dynamifyDSG<dynamic_symmetric_vertex, pbbs::empty, Graph>(GA); //Dynamify inside main
  HSet* h;

  if (type == 0) {
    h = new HSetDynArr(&_dynG);
    std::cout << "DYN_ARR VERSION\n" << std::endl;
  }
  else if (type == 1) {
    h = new HSetThreshold(&_dynG, GA.n);
    std::cout << "THRESHOLD VERSION\n" << std::endl;
  }

  clock.start();

  pbbs::sequence<pbbs::sequence<std::pair<uintE, uintE>>> batches = pbbs::sequence<pbbs::sequence<std::pair<uintE, uintE>>>(size);
  for (int i = 0; i < size; i++) {
    std::cout << "\n-----CHANGE " << (i + 1) << "-----" << std::endl;
    batches[i] = barabasi_albert::generate_updates(rand() % 1000 + 50, rand() % 100 + 5);
    h->insertEdges(getEdges(batches[i]));

    std::cout << "h-index: " << h->hindex << std::endl;
  }

  for (int i = size - 1; i >= 0; i--) {
    std::cout << "\n-----CHANGE " << (2 * size - i) << "-----" << std::endl;
    batches[i] = barabasi_albert::generate_updates(rand() % 1000 + 50, rand() % 100 + 5);
    h->eraseEdges(getEdges(batches[i]));

    std::cout << "h-index: " << h->hindex << std::endl;
  }

  std::cout << "\n\n------------------------------------" << endl;
  std::cout << "RUN TIME: " << clock.stop() << std::endl;

  return 0;
}


generate_symmetric_main(AppSubgraphCounting_runner, false);



























//-----DEBUGGING RUNNERS-----//

/* HSET DYN ARR
template <class Graph>
double AppHIndex_runner(Graph& GA, commandLine P) {
  std::cout << "### Application: H Index" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  //std::cout << "### Params: -k = " << k << " -e (epsilon) = " << epsilon << std::endl;
  std::cout << "### ------------------------------------" << endl;

  assert(P.getOption("-s"));
  // First, process the static graph -- construct the initial h-index structure
  // Then, assume we have some array that denotes batch edge insertions/deletions
  // In a serial for loop, process these batch updates, dynamically, including graph updates + h-index updates

  //EXAMPLE: ./h-index -s -rounds 1 "graph_test_3.txt"
  //         ./h-index -s -rounds 1 "graph_test_4.txt"

  //-rounds 1 to run once

  //------------------INITIALIZATION------------------//
  std::cout << "\n-----INIT GRAPH-----" << endl;
  auto _dynG = dynamifyDSG<dynamic_symmetric_vertex, pbbs::empty, Graph>(GA); //Dynamify inside main
  
  //HSet* h = new HSetDynArr(&_dynG);
  HSetDynArr* h = new HSetDynArr(&_dynG); //Temporarily using concrete type for debugging only

  //std::cout << "h-index: " << h->hindex << "    " << h->getH().size() << endl; //|H| = 2
  //std::cout << "B size: " << h->bSize << endl; //|B| = 1

  //Run Data for batches.size() = 10000 
  //real    2m46.280s
  //user    19m28.278s
  //sys	    1m34.892s
  pbbs::sequence<pbbs::sequence<std::pair<uintE, uintE>>> batches = pbbs::sequence<pbbs::sequence<std::pair<uintE, uintE>>>(100);

  for (int i = 0; i < batches.size(); i++) {
    //std::cout << "\n-----CHANGE " << (i + 1) << "-----" << endl;

    batches[i] = barabasi_albert::generate_updates(rand() % 1000 + 50, rand() % 100 + 5);

    h->insertEdges(batches[i]);
    //std::cout << "h-index: " << h->hindex << "    " << h->getH().size() << endl;
    //std::cout << "B size: " << h->bSize << endl;

    size_t cSize = 0;
    if (h->C.A[h->hindex] != nullptr) cSize = h->C.A[h->hindex]->size;

    //std::cout << "C[hindex] size: " << cSize << endl;

    size_t aboveH = 0;
    for (size_t j = 0; j < h->G->n; j++) {
      if (h->G->get_vertex(j).degree >= h->hindex) {
        aboveH++;
      }
    }
    //std::cout << "above H: " << aboveH << endl;
    assert(h->hindex + (cSize - h->bSize) == aboveH); //CHECK

  }

  for (int i = batches.size() - 1; i >= 0; i--) {
    //std::cout << "\n-----CHANGE " << (2 * batches.size() - i) << "-----" << endl;

    h->eraseEdges(batches[i]);

    //std::cout << "h-index: " << h->hindex << "    " << h->getH().size() << endl;
    //std::cout << "B size: " << h->bSize << endl;

    size_t cSize = 0;
    if (h->C.A[h->hindex] != nullptr) cSize = h->C.A[h->hindex]->size;

    //std::cout << "C[hindex] size: " << cSize << endl;

    size_t aboveH = 0;
    for (size_t j = 0; j < h->G->n; j++) {
      if (h->G->get_vertex(j).degree >= h->hindex) {
        aboveH++;
      }
    }
    //std::cout << "above H: " << aboveH << endl;

    assert(h->hindex + (cSize - h->bSize) == aboveH);
    
  }

  return 0;

}
*/

/*
// HSET HYBRID AKA HSET2
template <class Graph>
double AppHIndex_runner(Graph& GA, commandLine P) {
  std::cout << "### Application: H Index" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  //std::cout << "### Params: -k = " << k << " -e (epsilon) = " << epsilon << std::endl;
  std::cout << "### ------------------------------------" << endl;

  //assert(P.getOption("-s"));
  // First, process the static graph -- construct the initial h-index structure
  // Then, assume we have some array that denotes batch edge insertions/deletions
  // In a serial for loop, process these batch updates, dynamically, including graph updates + h-index updates

  //EXAMPLE: ./h-index -s -rounds 1 "graph_test_3.txt"
  //         ./h-index -s -rounds 1 "graph_test_4.txt"

  //-rounds 1 to run once

  //------------------INITIALIZATION------------------//
  std::cout << "\n-----INIT GRAPH-----" << endl;
  auto _dynG = dynamifyDSG<dynamic_symmetric_vertex, pbbs::empty, Graph>(GA); //Dynamify inside main
  
  //HSet* h = new HSetDynArr(&_dynG);
  HSetThreshold* h = new HSetThreshold(&_dynG, 100); //Temporarily using concrete type for debugging only

  //std::cout << "h-index: " << h->hindex << "    " << h->getH().size() << endl; //|H| = 2
  //std::cout << "B size: " << h->bSize << endl; //|B| = 1

  //Run data for batches.size() = 10,000 and threshold = 500
  //real	2m55.282s
  //user	20m38.886s
  //sys		1m36.433s

  //Run data for batches.size() = 10,000 and threshold = 100
  //real	3m9.755s
  //user	22m6.130s
  //sys		1m41.089s


  pbbs::sequence<pbbs::sequence<std::pair<uintE, uintE>>> batches = pbbs::sequence<pbbs::sequence<std::pair<uintE, uintE>>>(100);

  for (int i = 0; i < batches.size(); i++) {
    //std::cout << "\n-----CHANGE " << (i + 1) << "-----" << endl;

    batches[i] = barabasi_albert::generate_updates(rand() % 1000 + 50, rand() % 100 + 5);

    h->insertEdges(getEdges(batches[i]));

    uintE actualHSize = h->getH().size();
    assert(h->hindex == actualHSize);
    //std::cout << "h-index: " << h->hindex << "    " << actualHSize << endl;
    //std::cout << "B size: " << h->bSize << endl;

    size_t cSize = 0;
    if (h->getC(h->hindex) != nullptr) cSize = h->getC(h->hindex)->size;

    //std::cout << "C[hindex] size: " << cSize << endl;

    size_t aboveH = 0;
    for (size_t j = 0; j < h->G->n; j++) {
      if (h->G->get_vertex(j).degree >= h->hindex) {
        aboveH++;
      }
    }
    //std::cout << "above H: " << aboveH << endl;

    assert(h->hindex + (cSize - h->bSize) == aboveH);

  }

  for (int i = batches.size() - 1; i >= 0; i--) {
    std::cout << "\n-----CHANGE " << (2 * batches.size() - i) << "-----" << endl;

    h->eraseEdges(getEdges(batches[i]));

    uintE actualHSize = h->getH().size();
    assert(h->hindex == actualHSize);
    //std::cout << "h-index: " << h->hindex << "    " << actualHSize << endl;
    //std::cout << "B size: " << h->bSize << endl;

    size_t cSize = 0;
    if (h->getC(h->hindex) != nullptr) cSize = h->getC(h->hindex)->size;

    //std::cout << "C[hindex] size: " << cSize << endl;

    size_t aboveH = 0;
    for (size_t j = 0; j < h->G->n; j++) {
      if (h->G->get_vertex(j).degree >= h->hindex) {
        aboveH++;
      }
    }
    //std::cout << "above H: " << aboveH << endl;

    assert(h->hindex + (cSize - h->bSize) == aboveH);
    
  }

  return 0;
}
*/
