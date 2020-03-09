#include "h-index.h"

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
