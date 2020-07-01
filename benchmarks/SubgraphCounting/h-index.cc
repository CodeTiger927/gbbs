#include "h-index.h"
#include <iostream>

/*
int main() {

  symmetric_graph<symmetric_vertex, pbbs::empty> G = gbbs_io::read_unweighted_symmetric_graph("graph_test.txt", false);
  HSet h = HSet(G);
  h.threshold = 100;

  auto batch = pbbs::sequence<uintE>(6);
  par_for(0, 6, [&] (size_t i) {
    batch[i] = i;
  });

  h.erase(batch);

}
*/

/* Test rslice, interval from a to b in reverse (where a > b) = rslice(size - a, size - b)
int main() {
  auto s = pbbs::sequence<uintE>(5);

  for (int i = 0; i < 5; i++) {
    s[i] = i;
  }
  size_t a = 2;
  size_t b = 5;
  pbbs::sequence<uintE> x = s.rslice(a,b);

  cout << "A: " << a << endl;
  cout << "B: " << b << endl;

  for (int i = 0; i < x.size(); i++) {
    cout << x[i] << endl;
  }

}
*/
/*
int main() {

  // TODO: use integrated ligra main
  // e.g.,
  //template <class Graph>r
  //double AppKCore_runner(Graph& GA, commandLine P) {
  //}
  //generate_symmetric_main(AppKCore_runner, false);

  // HSet test
  cout << "HSet TEST" << endl;

  HSet h = HSet(1000, 100);
  *h.G = gbbs_io::read_unweighted_symmetric_graph("graph_test.txt", false);
  
  cout << h.insert(0) << endl;
  cout << h.insert(1) << endl;
  cout << h.insert(2) << endl;
  cout << h.insert(3) << endl;
  cout << h.insert(4) << endl;
  cout << h.insert(5) << endl << endl << endl;
  

  // graph test
  cout << "graph TEST" << endl;

  graph G;

  G.insertV(0);
  G.insertV(1);
  G.insertV(2);
  G.insertV(3);
  cout << G.connect(0,1) << endl;
  cout << G.connect(0,2) << endl;
  cout << G.connect(0,3) << endl;
  cout << G.connect(1,2) << endl;
  cout << G.connect(1,3) << endl;
  cout << G.connect(2,3) << endl;

  return 0;
}
*/


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

  //EXAMPLE: ./h-index -s -rounds 1 "graph_test.txt"

  //-rounds 1 to run once

  HSet<Graph> h = HSet<Graph>(GA);
  
  auto batch1 = pbbs::sequence<uintE>(1);
  batch1[0] = 4;

  auto batch2 = pbbs::sequence<uintE>(2);
  batch2[0] = 1;
  batch2[1] = 5;

  auto batch3 = pbbs::sequence<uintE>(1);
  batch3[0] = 0;

  auto batch4 = pbbs::sequence<uintE>(3);
  batch4[0] = 1;
  batch4[1] = 5;
  batch4[2] = 0;
  
  h.insert(batch1); //1 vertex with deg 1
  std::cout << "HINDEX: " << h.hindex << endl;
  std::cout << "B SIZE: " << h.B << endl;

  h.insert(batch3); //1 vertex with deg 3
  std::cout << "HINDEX: " << h.hindex << endl;
  std::cout << "B SIZE: " << h.B << endl;

  h.insert(batch2); //2 vertices with deg 2
  std::cout << "HINDEX: " << h.hindex << endl;
  std::cout << "B SIZE: " << h.B << endl;
  std::cout << "------------------------------------------" << endl;

  std::tuple<bool, sparse_table<uintE, pbbs::empty, hash_uintE>> info = h.containedInH(0);
  sparse_table<uintE, pbbs::empty, hash_uintE> H = std::get<1>(info);
  std::cout << "Contains Vertex 0: " << std::get<0>(info) << endl;
  pbbs::sequence<std::tuple<uintE, pbbs::empty>> flatten = H.entries();

  std::cout << "SIZE OF H: " << flatten.size() << endl;
  for (int i = 0; i < flatten.size(); i++) {
    cout << std::get<0>(flatten[i]) << endl;
  }

  cout << "------------------------------------------" << endl;
  h.erase(batch3); //removes 1 vertex with deg 3
  std::cout << "HINDEX: " << h.hindex << endl;
  std::cout << "B SIZE: " << h.B << endl;

  h.erase(batch2); //removes deg 2 and deg 3
  std::cout << "HINDEX: " << h.hindex << endl;
  std::cout << "B SIZE: " << h.B << endl;
  
  /*
  auto batch = pbbs::sequence<uintE>(6);
  par_for(0, 6, [&] (size_t i) {
    batch[i] = i;
  });

  h.insert(batch);
  std::cout << "HINDEX: " << h.hindex << endl;
  std::cout << "B SIZE: " << h.B << endl;
  std::cout << "------------------------------------------" << endl;

  std::tuple<bool, sparse_table<uintE, pbbs::empty, hash_uintE>> info = h.containedInH(4);
  sparse_table<uintE, pbbs::empty, hash_uintE> H = std::get<1>(info);
  std::cout << "Contains Vertex 0: " << std::get<0>(info) << endl;
  pbbs::sequence<std::tuple<uintE, pbbs::empty>> flatten = H.entries();

  std::cout << "SIZE OF H: " << flatten.size() << endl;
  for (int i = 0; i < flatten.size(); i++) {
    std::cout << std::get<0>(flatten[i]) << endl;
  }

  cout << "------------------------------------------" << endl;
  h.erase(batch);
  std::cout << "HINDEX: " << h.hindex << endl;
  std::cout << "B SIZE: " << h.B << endl;
  */
  return 0;
}

generate_symmetric_main(AppHIndex_runner, false);
