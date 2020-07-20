#include "h-index.h"
#include <iostream>
#include "ligra/pbbslib/dyn_arr.h"

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

/*
int main() {
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

  //EXAMPLE: ./h-index -s -rounds 1 "graph_test_3.txt"

  //-rounds 1 to run once

  //------------------INITIALIZATION------------------//
  std::cout << "\n-----INIT GRAPH-----" << endl;
  auto _dynG = dynamifyDSG<dynamic_symmetric_vertex, pbbs::empty, Graph>(GA); //Dynamify inside main
  HSet<Graph> h = HSet<Graph>(&_dynG);

  pbbs::sequence<uintE> vBatch = pbbs::sequence<uintE>(0); //General batch for vertices
  pbbs::sequence<std::pair<uintE, uintE>> eBatch = pbbs::sequence<std::pair<uintE, uintE>>(0); //General batch for edges

  std::cout << "h-index: " << h.hindex << endl; //|H| = 2
  std::cout << "B size: " << h.B.size() << endl; //|B| = 1



  //------------------CHANGE 1------------------//
  std::cout << "\n-----CHANGE 1-----" << endl;

  //Remove from HSet vertices that will be changed
  vBatch = pbbs::sequence<uintE>(2); //Temporarily removes vertices 1 and 3
  vBatch[0] = 1; vBatch[1] = 3;
  h.erase(vBatch); 

  //Add new vertices to graph
  vBatch = pbbs::sequence<uintE>(2); //Adds vertices 4 and 5 to graph
  vBatch[0] = 4; vBatch[1] = 5;
  h.G->batchAddVertices(vBatch);
  
  //Add new edges to graph
  eBatch = pbbs::sequence<std::pair<uintE, uintE>>(3); //Adds 1--3, 1--4, 3--5
  eBatch[0] = std::make_pair(1, 3); eBatch[1] = std::make_pair(1, 4); eBatch[2] = std::make_pair(3, 5);
  h.G->batchAddEdges(eBatch);

  //Add vertices to HSet
  vBatch = pbbs::sequence<uintE>(4); //Adds 4 and 5; Re-adds 1 and 3
  vBatch[0] = 1; vBatch[1] = 3; vBatch[2] = 4; vBatch[3] = 5; 
  h.insert(vBatch);

  std::cout << "h-index: " << h.hindex << endl; //|H| = 3
  std::cout << "B size: " << h.B.size() << endl; //|B| = 2


  //------------------CHANGE 2------------------//
  std::cout << "\n-----CHANGE 2-----" << endl;
  
  //Remove vertices from HSet
  vBatch = pbbs::sequence<uintE>(5); //Temporarily remove 0, 2, 4, 5; permanently removes 3
  vBatch[0] = 0; vBatch[1] = 2; vBatch[2] = 3; vBatch[3] = 4; vBatch[4] = 5;
  h.erase(vBatch);

  
  //Removes edges from graph
  eBatch = pbbs::sequence<std::pair<uintE, uintE>>(3); //Removes 1--3, 2--3, 3--5
  eBatch[0] = std::make_pair(1, 3); eBatch[1] = std::make_pair(2, 3); eBatch[2] = std::make_pair(3, 5);
  h.G->batchRemoveEdges(eBatch);

  //Remove vertices from graph
  vBatch = pbbs::sequence<uintE>(1); //Removes vertex 3
  vBatch[0] = 3;
  h.G->batchRemoveVertices(vBatch);


  //Add new edges to graph
  eBatch = pbbs::sequence<std::pair<uintE, uintE>>(3); //Adds edges 0--4, 2--4, 2--5, 4--5
  eBatch[0] = std::make_pair(2, 4); eBatch[1] = std::make_pair(2, 5); 
  eBatch[2] = std::make_pair(4, 5);
  h.G->batchAddEdges(eBatch);


  //Re-add changed vertices to HSet
  vBatch = pbbs::sequence<uintE>(4);
  vBatch[0] = 0; vBatch[1] = 2; vBatch[2] = 4; vBatch[3] = 5;
  h.insert(vBatch);

  std::cout << "h-index: " << h.hindex << endl; //|H| = 3
  std::cout << "B size: " << h.B.size() << endl; //|B| = 1

  //------------------CHANGE 3------------------//
  std::cout << "\n-----CHANGE 3-----" << endl;

  vBatch = pbbs::sequence<uintE>(5);
  vBatch[0] = 0; vBatch[1] = 1; vBatch[2] = 2; vBatch[3] = 4; vBatch[4] = 5;
  h.erase(vBatch);

  h.G->batchRemoveVertices(vBatch);

  std::cout << "h-index: " << h.hindex << endl; //|H| = 0
  std::cout << "B size: " << h.B.size() << endl; //|B| = 0

  return 0;
}

generate_symmetric_main(AppHIndex_runner, false);

