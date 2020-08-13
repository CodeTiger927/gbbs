#include "h-index-Alex.h"
#include "TriangleCounting.h"

#include <iostream>
#include <vector>
#include <fstream>


bool adj[1005][1005];

int main() {
  std::vector<std::pair<uintE,uintE>> edges;
  std::vector<std::pair<uintE,uintE>> edges2;

  std::ifstream in("graph.txt");
  for(int i = 0;i < 25571;++i) {
    uintE a,b;
    in >> a >> b;
    if(adj[a][b]) continue;
    adj[a][b] = true;
    adj[b][a] = true;
    if(a != b) {
      if(i < 25400) {
        edges.push_back(std::make_pair(a,b));
      }else{
        edges2.push_back(std::make_pair(a,b));
      }
    }
  }
  // cout << edges[0].first << "  " << edges[0].second << endl;

  sequence<std::pair<uintE,uintE>> allEdges = sequence<std::pair<uintE,uintE>>(edges.size(),[&](size_t i){return edges[i];});
  sequence<std::pair<uintE,uintE>> allEdges2 = sequence<std::pair<uintE,uintE>>(edges2.size(),[&](size_t i){return edges2[i];});

  // sequence<uintE> s = sequence<uintE>(1005,[&](size_t i){return i;});
  // cout << "hi" << endl;
  // dsg.batchAddVertices(s);
  // cout << "YAY" << endl;
  // dsg.batchAddEdges(allEdges);
  // cout << "YAY" << endl;
  // uintE maxD = 0;
  // for(int i = 0;i < 1005;++i) {
  //   if(dsg.v_data.A[i].degree == 546) cout << i << endl;
  //   maxD = std::max(maxD,dsg.v_data.A[i].degree);
  // }
  // cout << maxD << endl;
  std::cout << "Start Loading" << std::endl;
  initialize(1005);
  std::cout << "Done Loading" << std::endl;
  addEdges(allEdges);
  std::cout << total << std::endl;
  addEdges(allEdges2);
  std::cout << total << std::endl;
  removeEdges(allEdges2);
  std::cout << total << std::endl;
  removeEdges(allEdges);
  std::cout << total << std::endl;


  // std::cout << wedges[1][5] << endl;
  // uintE haha = wedges[1][5];
  // sequence<uintE> allH = hset.allH();
  // for(int i = 0;i < allH.size();++i) {
  //   if(dsg.existEdge(1,allH[i]) && dsg.existEdge(allH[i],5)) haha++;
  // }
  // std::cout << haha << std::endl;
}




// int main() {
//   std::vector<std::pair<uintE,uintE>> v1 = {std::make_pair(0,1),std::make_pair(0,2),std::make_pair(1,2)};
//   std::vector<std::pair<uintE,uintE>> v2 = {std::make_pair(3,4),std::make_pair(3,5),std::make_pair(3,6),std::make_pair(4,5),std::make_pair(4,6),std::make_pair(5,6)};
//   std::vector<std::pair<uintE,uintE>> v3 = {std::make_pair(0,3),std::make_pair(0,4),std::make_pair(0,5),std::make_pair(0,6),std::make_pair(1,3),std::make_pair(1,4),std::make_pair(1,5),std::make_pair(1,6),std::make_pair(2,3),std::make_pair(2,4),std::make_pair(2,5),std::make_pair(2,6)};

//   std::vector<std::pair<uintE,uintE>> v4 = {std::make_pair(0,1)};
//   std::vector<std::pair<uintE,uintE>> v5 = {std::make_pair(0,2),std::make_pair(0,3),std::make_pair(1,2),std::make_pair(1,3)};
//   std::vector<std::pair<uintE,uintE>> v6 = {std::make_pair(2,3)};

//   sequence<std::pair<uintE,uintE>> s = sequence<std::pair<uintE,uintE>>(3,[&](size_t i){return v1[i];});
//   sequence<std::pair<uintE,uintE>> s2 = sequence<std::pair<uintE,uintE>>(6,[&](size_t i){return v2[i];});
//   sequence<std::pair<uintE,uintE>> s10 = sequence<std::pair<uintE,uintE>>(12,[&](size_t i){return v3[i];});

//   sequence<std::pair<uintE,uintE>> s12 = sequence<std::pair<uintE,uintE>>(1,[&](size_t i){return v4[i];});
//   sequence<std::pair<uintE,uintE>> s13 = sequence<std::pair<uintE,uintE>>(4,[&](size_t i){return v5[i];});
//   sequence<std::pair<uintE,uintE>> s14 = sequence<std::pair<uintE,uintE>>(1,[&](size_t i){return v6[i];});

//   sequence<std::pair<uintE,uintE>> s5 = sequence<std::pair<uintE,uintE>>(2,[&](size_t i){return std::make_pair(i,2);});


//   // sequence<std::pair<uintE,uintE>> s2 = sequence<std::pair<uintE,uintE>>(5,[&](size_t i){return std::std::make_pair(i,6);});
//   sequence<uintE> s3 = sequence<uintE>(10,[&](size_t i){return i;});
//   sequence<std::pair<uintE,uintE>> s4 = sequence<std::pair<uintE,uintE>>(10,[&](size_t i){return std::make_pair(i,0);});
//   sequence<uintE> s6 = sequence<uintE>(3,[&](size_t i){return i;});
//   sequence<uintE> s9 = sequence<uintE>(4,[&](size_t i){return i + 3;});
//   sequence<std::pair<uintE,uintE>> s8 = sequence<std::pair<uintE,uintE>>(3,[&](size_t i){return std::make_pair(i,2);});
//   sequence<std::pair<uintE,uintE>> s7 = sequence<std::pair<uintE,uintE>>(4,[&](size_t i){return std::make_pair(i + 3,3);});
//   sequence<std::pair<uintE,uintE>> s11 = sequence<std::pair<uintE,uintE>>(7,[&](size_t i){return std::make_pair(i,6);});

//   hset.resizeV(10);
//   hset.insert(s4);

//   dsg.batchAddVertices(s3);


//   // hset.modify(s8);
//   // hset.modify(s7);
//   // cout << hset.BSize << " " << hset.HIndex << endl;
//   // hset.modify(s11);
//   // cout << hset.BSize << " " << hset.HIndex << endl;


//   // // cout << "haha " << total << endl;
//   addEdges(s);
//   cout << total << endl;
//   addEdges(s2);
//   cout << total << endl;
//   addEdges(s10);
//   cout << total << endl;

//   removeEdges(s10);
//   cout << total << endl;
//   removeEdges(s2);
//   cout << total << endl;
//   removeEdges(s);
//   cout << total << endl;
//   // cout << total << endl;
//   // cout << hset.HIndex << " " << hset.BSize << endl;
//   // cout << hset.allH()[0] << " " << hset.allH()[1] << endl;
//   // cout << "LOL WHAT "<< hset.HIndex << endl;

//   // cout << wedges[0][1] << endl;

//   // sequence<uintE> hnodes = hset.allH();
//   // for(int i = 0;i < hnodes.size();++i) {
//   //   cout << hnodes[i] << endl;
//   // }
//   // cout << hnodes.size() << endl;
//   // cout << hset.HIndex << endl;
//   // removeEdges(s2);
//   // cout << total << endl;
//   // removeEdges(s10);
//   // cout << total << endl;
//   // removeEdges(s);
//   // cout << total << endl;

//   // sequence<uintE> hehe = hset.allH();
//   // cout << hehe[0] << "   " << hehe[1] << endl;
//   // cout << hset.HIndex << endl;
//   // removeEdges(s);
//   // cout << total << endl;







//   // cout << (triObj(1,6,3) == triObj(3,6,1)) << endl;


//   // cout << hs.B.find(0,false) << endl;
// }





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

// template <class Graph>
// double AppHIndex_runner(Graph& GA, commandLine P) {
//   std::cout << "### Application: H Index" << std::endl;
//   std::cout << "### Graph: " << P.getArgument(0) << std::endl;
//   std::cout << "### Threads: " << num_workers() << std::endl;
//   std::cout << "### n: " << GA.n << std::endl;
//   std::cout << "### m: " << GA.m << std::endl;
//   //std::cout << "### Params: -k = " << k << " -e (epsilon) = " << epsilon << std::endl;
//   std::cout << "### ------------------------------------" << endl;

//   assert(P.getOption("-s"));
//   // First, process the static graph -- construct the initial h-index structure
//   // Then, assume we have some array that denotes batch edge insertions/deletions
//   // In a serial for loop, process these batch updates, dynamically, including graph updates + h-index updates

//   //EXAMPLE: ./h-index -s -rounds 1 "graph_test.txt"

//   //-rounds 1 to run once

//   //For some reason GA has type symmetric_graph<csv_bytepd_amortized, pbbs::empty> which doesn't match my symmetric_graph<symmetric_vertex, pbbs::empty>&
//   //I don't know how to fix that so for testing I just commented out 73 and uncommented out line 69-70
//   // TODO: This part needs to be fixed -- why are you re-reading the graph? You have to
//   // templatize everything like we discussed.
//   symmetric_graph<symmetric_vertex, pbbs::empty> G = gbbs_io::read_unweighted_symmetric_graph(P.getArgument(0), false);
//   HSet<Graph> h = HSet(GA);
//   h.threshold = 1; //For debugging, to make sure that whatever threshold is, it still works
  
//   //HSet h = HSet(GA); 
//   auto batch = pbbs::sequence<uintE>(6);
//   par_for(0, 6, [&] (size_t i) {
//     batch[i] = i;
//   });
  
//   h.insert(batch);
//   cout << "HINDEX: " << h.hindex << endl;

//   return 0;
// }

// generate_symmetric_main(AppHIndex_runner, false);
