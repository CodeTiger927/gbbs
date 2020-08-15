// TODO from myself:
// 1. Clean up everything. Wayyyy too many repetitive code, and compiling takes so long right now
// 2. A lot of the initializations are not even automatic, and need to be manually changed
//   a. The initialization function doesn't always work
//   b. For some of the header files, you need to manually change the Size Graph first, which is not ideal
//   c. I don't even know why there are 2 variables controlling the size of the graph in h-index-Alex.h
//   There's probably more, but this is what i remember for now.
// 3. Similarly, I needed to manually figure out what the max degree is, or else the allH function will probably time out horribly.
// 4. Speaking of which, we need a more efficient allH function that is not O(N) worst case.
// 5. Also more testing in general? There probably are a lot more bugs that I missed.
// 6. I currently uses a lot of unnecessary memory that can probably be much further reduced(I wouldn't be surprised if we cut down 5 times the memory)
// 7. Similarly, I used a lot of algorithms that probably aren't optimal for performance, so we can change those too.

#include "h-index-Alex.h"
#include "TriangleCounting.h"

#include <iostream>
#include <vector>
#include <fstream>


bool adj[1005][1005];

int main() {
  std::vector<std::pair<uintE,uintE>> edges;
  std::vector<std::pair<uintE,uintE>> edges2;

  std::vector<std::pair<uintE,uintE>> edges3;
  std::vector<std::pair<uintE,uintE>> edges4;


  std::ifstream in("graph.txt");
  for(int i = 0;i < 25571;++i) {
    uintE a,b;
    in >> a >> b;
    if(adj[a][b]) continue;
    adj[a][b] = true;
    adj[b][a] = true;
    if(a != b) {
      if(i < 15400) {
        edges.push_back(std::make_pair(a,b));
      }else{
        edges2.push_back(std::make_pair(a,b));
      }

      if(i % 2) {
        edges3.push_back(std::make_pair(a,b));
      }else{
        edges4.push_back(std::make_pair(a,b));
      }
    }
  }

  sequence<std::pair<uintE,uintE>> allEdges = sequence<std::pair<uintE,uintE>>(edges.size(),[&](size_t i){return edges[i];});
  sequence<std::pair<uintE,uintE>> allEdges2 = sequence<std::pair<uintE,uintE>>(edges2.size(),[&](size_t i){return edges2[i];});
  
  sequence<std::pair<uintE,uintE>> allEdges3 = sequence<std::pair<uintE,uintE>>(edges3.size(),[&](size_t i){return edges3[i];});
  sequence<std::pair<uintE,uintE>> allEdges4 = sequence<std::pair<uintE,uintE>>(edges4.size(),[&](size_t i){return edges4[i];});

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
  removeEdges(allEdges3);
  std::cout << total << std::endl;
  removeEdges(allEdges4);
  std::cout << total << std::endl;


  // std::cout << wedges[1][5] << endl;
  // uintE haha = wedges[1][5];
  // sequence<uintE> allH = hset.allH();
  // for(int i = 0;i < allH.size();++i) {
  //   if(dsg.existEdge(1,allH[i]) && dsg.existEdge(allH[i],5)) haha++;
  // }
  // std::cout << haha << std::endl;
}


