#include <cstring>
#include <stdlib.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <fstream>
#include <cmath>

#include "pbbslib/sequence_ops.h"

#include "pbbslib/merge_sort.h"

#include "pbbslib/monoid.h"

#include "ligra/pbbslib/sparse_table.h"

#include "ligra/pbbslib/dyn_arr.h"


size_t INIT_DYN_GRAPH_EDGE_SIZE = 4;

std::function<void()> get_deletion_fn(void*, void*);
std::function<void()> get_deletion_fn(void*, void*, void*);
std::function<void()> get_deletion_fn(void*, void*, void*, void*);



struct hash_uintE {
  inline size_t operator () (const uintE & a) {return pbbs::hash64_2(a);}
};

struct dynamic_vertex_data {
  uintE degree;  // vertex degree
  sparse_table<uintE,bool,hash_uintE> neighbors; // Sparse_table
  uintE stored; // Sparse_table
};

template <class W>
struct dynamic_symmetric_vertex {
  using vertex = dynamic_symmetric_vertex<W>;
  using edge_type = std::tuple<uintE, W>;
  using edges_type = sparse_table<uintE,bool,hash_uintE>;

  uintE degree;
  uintE stored;
  edges_type* neighbors = nullptr;

  dynamic_symmetric_vertex(dynamic_vertex_data& vdata) {
    neighbors = vdata.neighbors;
    degree = vdata.degree;
    stored = vdata.stored;
  }

  // Return all the neighbors of that vertex
  sequence<uintE> getAllNeighbor() {
    auto entries = neighbors -> entries();
    sequence<uintE> res = sequence<uintE>(entries.size(),[&](size_t i) {
      std::tuple<uintE,bool> cur = entries[i];
      return std::get<0>(cur);
    });
    return res;
  }

  // Clear
  void clear() { neighbors -> del(); }
};

template <template <class W> class vertex_type, class W>
struct dynamic_symmetric_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  using edge_type = typename vertex::edge_type;

  pbbslib::dyn_arr<dynamic_vertex_data> v_data;

  /* number of vertices in G */
  size_t n;
  /* number of edges in G */
  size_t m;

  // Boolean array that keeps track if a vertex exists
  pbbslib::dyn_arr<bool> existVertices;

  /* called to delete the graph */
  std::function<void()> deletion_fn;

  dynamic_symmetric_graph(pbbslib::dyn_arr<dynamic_vertex_data> & _v_data, size_t _n, size_t _m,
    std::function<void()> _deletion_fn,pbbslib::dyn_arr<bool> & _existVertices)
      : n(_n),
        m(_m),
        deletion_fn(_deletion_fn) {
          v_data = pbbslib::dyn_arr<dynamic_vertex_data>(_v_data.size);
          par_for(0,_v_data.size,1,[&](size_t i) {
            v_data.A[i] = _v_data.A[i];
          });
          // existVertices = pbbslib::dyn_arr<bool>(_existVertices.size);
          // par_for(0,_existVertices.size,1,[&](size_t i) {
          //   existVertices.A[i] = _existVertices.A[i];
          // });
  }

  void del() {
    deletion_fn();
  }

  vertex get_vertex(uintE i) {
    return vertex(v_data.A[i]);
  }


  // A function to check whether or not it is worth resizing
  void adjustVdata(size_t amount) {
    amount = std::max(amount,INIT_DYN_GRAPH_EDGE_SIZE);
    size_t cur = v_data.capacity;
    if(amount != 0 && ceil(log2(amount)) == ceil(log2(cur))) {
      return;
    }
    size_t nC = 1 << ((size_t)(ceil(log2(amount))));
    dynamic_vertex_data* nA = pbbslib::new_array_no_init<dynamic_vertex_data>(nC);
    par_for(0,v_data.size,1,[&](size_t i){nA[i] = v_data.A[i];});
    pbbslib::free_array(v_data.A);
    v_data.A = nA;
    v_data.capacity = nC;

    bool* nAA = pbbslib::new_array_no_init<bool>(nC);
    par_for(0,existVertices.size,1,[&](size_t i){nAA[i] = existVertices.A[i];});
    pbbslib::free_array(existVertices.A);
    existVertices.A = nAA;
    existVertices.capacity = nC;
  }



  // A function that determines whether it is worth it to resize neighbors
  void adjustNeighbors(uintE u) {
    size_t amount = std::max((size_t)(v_data.A[u].stored),INIT_DYN_GRAPH_EDGE_SIZE);
    size_t cur = v_data.A[u].neighbors.m;
    if((amount << 1) <= cur && std::max(INIT_DYN_GRAPH_EDGE_SIZE,(size_t)v_data.A[u].degree) << 2 >= cur) return; 
    auto entries = v_data.A[u].neighbors.entries();
    // Sparsity of always at least 2, so find takes on average 2 steps to find empty
    v_data.A[u].neighbors = make_sparse_table<uintE,bool,hash_uintE>(2 * std::max(INIT_DYN_GRAPH_EDGE_SIZE,(size_t)v_data.A[u].degree),std::make_tuple(UINT_E_MAX,false),hash_uintE());
    par_for(0,entries.size(),1,[&](size_t i) {
      uintE cur = std::get<0>(entries[i]);
      if(std::get<1>(entries[i])) v_data.A[u].neighbors.insert(std::make_tuple(cur,true));
    });
    v_data.A[u].stored = v_data.A[u].degree;
  }

  // Add and initialize vertices in parallel
  // All indices in vertices need to be unique.
  void batchAddVertices(sequence<uintE> & vertices) {
    size_t size = vertices.size();
    uintE ma = pbbs::reduce(vertices,pbbs::maxm<uintE>());
    adjustVdata(std::max((size_t)ma,v_data.capacity));

    pbbs::sequence<uintE> existingVertices = pbbs::filter(vertices,[&](uintE i){return i >= this -> n;});
    // cout << existingVertices.size() << endl;
    size_t sumV = existingVertices.size();
    par_for(0,sumV,1,[&](size_t i) {
      uintE id = existingVertices[i];
      existVertices.A[id] = true;

      v_data.A[id].neighbors = make_sparse_table<uintE,bool,hash_uintE>(INIT_DYN_GRAPH_EDGE_SIZE,std::make_tuple(UINT_E_MAX,false),hash_uintE());
      v_data.A[id].degree = 0;
      v_data.A[id].stored = 0;
    });
    this -> n = std::max(this -> n,(size_t)ma);

  }




  // Check if an edge exists
  bool existEdge(uintE v,uintE u) {
    return v_data.A[v].neighbors.find(u,false);
  }

  // Add a series of edges in parallel
  // All edges need to be unique
  void batchAddEdges(uintE u,sequence<uintE> edges) {
    pbbs::sequence<uintE> ds = pbbs::filter(edges,[&](uintE i){return !existEdge(u,i);});
    uintE sumW = ds.size();

    this -> m += sumW;

    v_data.A[u].degree += ds.size();
    v_data.A[u].stored += ds.size();
    adjustNeighbors(u);

    par_for(0,ds.size(),1,[&](size_t i) {
      uintE v = ds[i];
      v_data.A[v].degree++;
      v_data.A[v].stored++;
      adjustNeighbors(v);

      v_data.A[v].neighbors.insert(std::make_tuple(u,true));
      v_data.A[u].neighbors.insert(std::make_tuple(v,true));
    });
  }


  // Add edges in the forms of pairs in parallel
  void batchAddEdges(sequence<std::pair<uintE,uintE>> edges) {
    sequence<std::pair<uintE,uintE>> ds = filter(edges,[&](std::pair<uintE,uintE> i){return !existEdge(i.first,i.second);});
    
    uintE sumW = ds.size();

    this -> m += sumW;

    sequence<std::pair<uintE,uintE>> allS = sequence<std::pair<uintE,uintE>>(2 * ds.size(),[&](size_t i) {
      if(i % 2) {
        return std::make_pair(ds[i / 2].first,ds[i / 2].second);
      }else{
        return std::make_pair(ds[i / 2].second,ds[i / 2].first);
      }
    });


    // Sorting to recalculate the degree

    // Using merge sort for stability, as it guarentees NlogN work regardless of data specifics.
    sequence<std::pair<uintE,uintE>> all = merge_sort(allS,[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a < b;});

    par_for(0,all.size(),[&](size_t i) {
      std::pair<uintE,uintE> cur = all[i];
      if(i != 0) {
        if(cur.first != all[i - 1].first) {
          v_data.A[all[i - 1].first].degree += i;
          v_data.A[all[i - 1].first].stored += i;
        }
      }
    });
    // cout << start.m << endl;
    v_data.A[all[2 * ds.size() - 1].first].degree += 2 * ds.size();
    v_data.A[all[2 * ds.size() - 1].first].stored += 2 * ds.size();

    par_for(0,all.size(),[&](size_t i) {
      std::pair<uintE,uintE> cur = all[i];
      if(i != 0) {
        if(cur.first != all[i - 1].first) {
          v_data.A[cur.first].degree -= i;
          v_data.A[cur.first].stored -= i;
          adjustNeighbors(cur.first);
        }
      }
    });
    
    adjustNeighbors(all[0].first);

    par_for(0,all.size(),[&](size_t i) {
      std::pair<uintE,uintE> cur = all[i];
      v_data.A[cur.first].neighbors.insert(std::make_tuple(cur.second,true));
    });

  }

  // Remove multiple vertices in parallel
  void batchRemoveVertices(sequence<uintE> & vertices) {
    size_t size = vertices.size();

    pbbs::sequence<uintE> existingVertices = pbbs::filter(vertices,[&](uintE i){return existVertices.A[i];});
    size_t sumV = existingVertices.size();
    par_for(0,sumV,1,[&](size_t i) {

      uintE id = existingVertices[i];
      existVertices.A[id] = false;

      v_data.A[id].neighbors.del();
      v_data.A[id].degree = 0;
    });

    adjustVdata(v_data.size - sumV);

    this -> n -= sumV;
  }

  void batchRemoveEdges(uintE u,sequence<uintE> edges) {
    pbbs::sequence<uintE> ds = pbbs::filter(edges,[&](uintE i){return existEdge(u,i);});
    uintE sumW = ds.size();

    this -> m -= sumW;

    v_data.A[u].degree -= ds.size();

    par_for(0,ds.size(),1,[&](size_t i) {
      uintE v = ds[i];

      // Remove u from v
      v_data.A[v].degree--;
      v_data.A[v].neighbors.change(u,false);
      v_data.A[u].neighbors.change(v,false);

      adjustNeighbors(v);
    });
    adjustNeighbors(u);
  }

    // Remove edges in the forms of pairs in parallel
  void batchRemoveEdges(sequence<std::pair<uintE,uintE>> edges) {
    sequence<std::pair<uintE,uintE>> ds = filter(edges,[&](std::pair<uintE,uintE> i){return existEdge(i.first,i.second);});
    uintE sumW = ds.size();

    this -> m -= sumW;
    sequence<std::pair<uintE,uintE>> allS = sequence<std::pair<uintE,uintE>>(2 * ds.size(),[&](size_t i) {
      if(i % 2) {
        return std::make_pair(ds[i / 2].first,ds[i / 2].second);
      }else{
        return std::make_pair(ds[i / 2].second,ds[i / 2].first);
      }
    });

    // Sorting to recalculate the degree

    // Using merge sort for stability, as it guarentees NlogN work regardless of data specifics.
    sequence<std::pair<uintE,uintE>> all = merge_sort(allS,[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a < b;});

    par_for(0,all.size(),[&](size_t i) {
      std::pair<uintE,uintE> cur = all[i];
      v_data.A[cur.first].neighbors.change(cur.second,false);
      if(i != 0) {
        if(cur.first != all[i - 1].first) {
          v_data.A[all[i - 1].first].degree -= i;
        }
      }
    });

    v_data.A[all[2 * ds.size() - 1].first].degree -= 2 * ds.size();

    adjustNeighbors(all[0].first);

    par_for(0,all.size(),1000000,[&](size_t i) {
      std::pair<uintE,uintE> cur = all[i];
      if(i != 0) {
        if(cur.first != all[i - 1].first) {
          v_data.A[cur.first].degree += i;
          adjustNeighbors(cur.first);
        }
      }
    });

  }
};

// Creates an empty dynamic symmetric graph
template <template <class W> class vertex_type, class W>
dynamic_symmetric_graph<dynamic_symmetric_vertex,W> createEmptyDynamicSymmetricGraph() {

  std::function<void()> doNothing = []() {};
  pbbslib::dyn_arr<dynamic_vertex_data> _vData = pbbslib::dyn_arr<dynamic_vertex_data>(1);
  pbbslib::dyn_arr<bool> tmp = pbbslib::dyn_arr<bool>(1);
  return dynamic_symmetric_graph<dynamic_symmetric_vertex,W>(_vData,0,0,doNothing,tmp);
}

// Creates a dynamic symmetric graph based on another graph
template <template <class W> class vertex_type, class W>
dynamic_symmetric_graph<dynamic_symmetric_vertex,W> dynamifyDSG(symmetric_graph<symmetric_vertex,pbbs::empty> & G) {

  dynamic_symmetric_graph<dynamic_symmetric_vertex,W> dsg = createEmptyDynamicSymmetricGraph<dynamic_symmetric_vertex,W>();

  pbbs::sequence<uintE> v = pbbs::sequence<uintE>(G.n,[&](size_t i){return i;});
  dsg.batchAddVertices(v);

  for(int i = 0;i < G.n;i++) {
    pbbs::sequence<uintE> s = pbbs::sequence<uintE>(G.get_vertex(i).degree,[&](size_t j){return G.get_vertex(i).getOutNeighbor(j);});
    dsg.batchAddEdges(i,s);
  }
  
  return dsg;
}