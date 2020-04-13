#include <iostream>
#include <map>
#include <set>
#include <limits>
#include <vector>
#include <math.h>
#include "ligra/ligra.h"
#include "dynamic_symmetric_graph.h"
#include "dyn_arr_sym_graph.h"
#include "ligra/pbbslib/sparse_table.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "pbbslib/sequence.h"

uintE N;

dynamic_symmetric_graph<dynamic_symmetric_vertex,uintE> dsg;
HSet hset = nullptr;
uintE hedges[N][N];
uintE total = 0;

bool checkTF(sparse_table<uintE,bool,hash_uintE> st,uintE v) {
  return st.find(v,false);
}


// Basic concept and not parallel
uintE connectEdge(uintE v,uintE u) {
  if(dsg.existEdge(u,v)) return total;
  hset.connect(u,v);

  auto entries = hset.P.entries();
  total += hedges[v][u];
  par_for(0,entries.size(),1,[&](size_t i){
    if(dsg.existEdge(std::get<0>(entries[0]),u) && dsg.existEdge(std::get<0>(entries[0]),v)) ++total;
  });
  if(!checkTF(hset.P,v)) {
    par_for(0,dsg.v_data[v].entries.size,1,[&](size_t i) {
      ++hedges[u][std::get<0>(dsg.v_data[v].entries[i])];
      ++hedges[std::get<0>(dsg.v_data[v].entries[i])][u];
    });
  }

  if(!checkTF(hset.P,u)) {
    par_for(0,dsg.v_data[u].entries.size,1,[&](size_t i) {
      ++hedges[v][std::get<0>(dsg.v_data[u].entries[i])];
      ++hedges[std::get<0>(dsg.v_data[u].entries[i])][v];
    });
  }
  dsg.addEdge(u,v);
  return total;
}