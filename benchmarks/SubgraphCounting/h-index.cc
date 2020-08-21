//#include "TriangleCounting.h"
#include "h-index.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "utils/generators/barabasi_albert.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <assert.h>


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
  HSetDynArr* h = new HSetDynArr(&_dynG); //Temporarily using concrete type for debugging only

  //std::cout << "h-index: " << h->hindex << "    " << h->getH().size() << endl; //|H| = 2
  //std::cout << "B size: " << h->bSize << endl; //|B| = 1

  //Run Data for batches.size() = 10000 
  //real    2m46.280s
  //user    19m28.278s
  //sys	    1m34.892s
  pbbs::sequence<pbbs::sequence<std::pair<uintE, uintE>>> batches = pbbs::sequence<pbbs::sequence<std::pair<uintE, uintE>>>(10000);

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

/*
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
  HSet2* h = new HSet2(&_dynG, 20); //Temporarily using concrete type for debugging only

  std::cout << "h-index: " << h->hindex << "    " << h->getH().size() << endl; //|H| = 2
  std::cout << "B size: " << h->bSize << endl; //|B| = 1

  pbbs::sequence<pbbs::sequence<std::pair<uintE, uintE>>> batches = pbbs::sequence<pbbs::sequence<std::pair<uintE, uintE>>>(7);

  for (int i = 0; i < batches.size(); i++) {
    std::cout << "\n-----CHANGE " << (i + 1) << "-----" << endl;

    batches[i] = barabasi_albert::generate_updates(rand() % 1000 + 50, rand() % 100 + 5);

    h->insertEdges(batches[i]);
    std::cout << "h-index: " << h->hindex << "    " << h->getH().size() << endl;
    std::cout << "B size: " << h->bSize << endl;

    size_t cSize = 0;
    if (h->getC(h->hindex) != nullptr) cSize = h->getC(h->hindex)->size;

    std::cout << "C[hindex] size: " << cSize << endl;

    size_t aboveH = 0;
    for (size_t j = 0; j < h->G->n; j++) {
      if (h->G->get_vertex(j).degree >= h->hindex) {
        aboveH++;
      }
    }
    cout << "above H: " << aboveH << endl;

    if (h->hindex + (cSize - h->bSize) != aboveH) cout << "-------------------VERY BAD------------------" << endl;

  }

  for (int i = batches.size() - 1; i >= 0; i--) {
    std::cout << "\n-----CHANGE " << (2 * batches.size() - i) << "-----" << endl;

    h->eraseEdges(batches[i]);

    std::cout << "h-index: " << h->hindex << "    " << h->getH().size() << endl;
    std::cout << "B size: " << h->bSize << endl;

    size_t cSize = 0;
    if (h->getC(h->hindex) != nullptr) cSize = h->getC(h->hindex)->size;

    std::cout << "C[hindex] size: " << cSize << endl;

    size_t aboveH = 0;
    for (size_t j = 0; j < h->G->n; j++) {
      if (h->G->get_vertex(j).degree >= h->hindex) {
        aboveH++;
      }
    }
    cout << "above H: " << aboveH << endl;

    if (h->hindex + (cSize - h->bSize) != aboveH) cout << "-------------------VERY BAD------------------" << endl;
    
  }


  return 0;
}
*/

generate_symmetric_main(AppHIndex_runner, false);
