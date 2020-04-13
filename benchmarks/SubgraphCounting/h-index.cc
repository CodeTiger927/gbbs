#include "h-index.h"
#include <iostream>




int main() {

  auto G = gbbs_io::read_unweighted_symmetric_graph("graph_test.txt", false);
  HSet h = HSet(2, G);
  sequence<uintE> batch = sequence<uintE>(6);
  par_for(0, 6, [&] (size_t i) {
    batch[i] = (2 * i) % 6;
  });

  cout << "New |H|" << h.insert(batch) << endl;
}

/*
int main() {

  // TODO: use integrated ligra main
  // e.g.,
  //template <class Graph>
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

  assert(P.getOption("-s"));
  // First, process the static graph -- construct the initial h-index structure
  // Then, assume we have some array that denotes batch edge insertions/deletions
  // In a serial for loop, process these batch updates, dynamically, including graph updates + h-index updates

  //Temporary Test for H-index
  //EXAMPLE: ./h-index -s "graph_test.txt"
  //Why does it run 3 times?
  // TODO: You have to pass -rounds 1 to make it run only 1 time (or however many times you would like).
  // The reason why it runs 3 times is because sometimes, the cache is cold on the first run;
  // the second run warms up the cache, and may give faster times. Usually we report the 
  // minimum of all 3 runs for timing tests.
  // TODO: Why is G being read in this way? GA is the passed in graph.
  symmetric_graph<symmetric_vertex, pbbs::empty> G = gbbs_io::read_unweighted_symmetric_graph(P.getArgument(0), false);
  // TODO: Threshold should be determined based on G, not written as a constant
  HSet h = HSet(100, G);
  
  cout << h.insert(0) << endl;
  cout << h.insert(1) << endl;
  cout << h.insert(2) << endl;
  cout << h.insert(3) << endl;
  cout << h.insert(4) << endl;
  cout << h.insert(5) << endl;
  cout << h.erase(0) << endl;
  cout << h.erase(1) << endl;
  cout << h.erase(2) << endl << endl << endl;

  cout << "Execution Complete" << endl;

  return 0;
}

generate_symmetric_main(AppHIndex_runner, false);
*/
