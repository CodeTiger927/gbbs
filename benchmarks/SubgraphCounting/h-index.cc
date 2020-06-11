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

  //For some reason GA has type symmetric_graph<csv_bytepd_amortized, pbbs::empty> which doesn't match my symmetric_graph<symmetric_vertex, pbbs::empty>&
  //I don't know how to fix that so for testing I just commented out 73 and uncommented out line 69-70
  symmetric_graph<symmetric_vertex, pbbs::empty> G = gbbs_io::read_unweighted_symmetric_graph(P.getArgument(0), false);
  HSet<Graph> h = HSet(GA);
  h.threshold = 1; //For debugging, to make sure that whatever threshold is, it still works
  
  //HSet h = HSet(GA); 
  auto batch = pbbs::sequence<uintE>(6);
  par_for(0, 6, [&] (size_t i) {
    batch[i] = i;
  });
  
  h.insert(batch);
  cout << "HINDEX: " << h.hindex << endl;

  return 0;
}

generate_symmetric_main(AppHIndex_runner, false);
