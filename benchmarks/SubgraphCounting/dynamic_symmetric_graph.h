#include <cstring>
#include <stdlib.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <fstream>
#include <cmath>

#include "pbbslib/sequence_ops.h"

#include "pbbslib/monoid.h"

#include "ligra/pbbslib/sparse_table.h"

#include "ligra/pbbslib/dyn_arr.h"

// Initial size for sparse_table
uintE INIT_VALUE_FOR_SIZE = 10;

std::function<void()> get_deletion_fn(void*, void*);
std::function<void()> get_deletion_fn(void*, void*, void*);
std::function<void()> get_deletion_fn(void*, void*, void*, void*);

struct hash_uintE {
  inline size_t operator () (const uintE & a) {return pbbs::hash64_2(a);}
};

struct dynamic_vertex_data {
  uintE degree;  // vertex degree
  sparse_table<uintE,bool,hash_uintE> neighbors; // Sparse_table
};

template <class W>
struct dynamic_symmetric_vertex {
  using vertex = dynamic_symmetric_vertex<W>;
  using edge_type = std::tuple<uintE, W>;
  using edges_type = sparse_table<uintE,bool,hash_uintE>;

  uintE degree;
  edges_type neighbors = make_sparse_table<uintE,bool,hash_uintE>(INIT_VALUE_FOR_SIZE,std::make_tuple(UINT_E_MAX,false),hash_uintE());

  dynamic_symmetric_vertex(dynamic_vertex_data& vdata) {
    neighbors = vdata.neighbors;
    degree = vdata.degree;
  }

  sequence<uintE> allNeighbor() {
    auto entries = this.neighbors.entries();
    sequence<uintE> res = sequence<uintE>(entries.size(),[&](size_t i) {
      std::tuple<uintE,bool> cur = entries[i];
      return std::get<0>(cur);
    });
    return res;
  }

  void clear() { neighbors.del(); }
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

  // Sparse table that keeps track if a vertex exists
  sparse_table<uintE,bool,hash_uintE> existVertices;

  /* called to delete the graph */
  std::function<void()> deletion_fn;

  dynamic_symmetric_graph(pbbslib::dyn_arr<dynamic_vertex_data> &_v_data, size_t n, size_t m,
    std::function<void()> _deletion_fn,sparse_table<uintE,bool,hash_uintE> _existVertices)
      : n(n),
        m(m),
        existVertices(_existVertices),
        deletion_fn(_deletion_fn) {
          this -> v_data = pbbslib::dyn_arr<dynamic_vertex_data>(_v_data.size);
          par_for(0,_v_data.size,1,[&](size_t i) {
            this -> v_data.A[i] = _v_data.A[i];
          });
  }

  // Find union of edges between two vertices
  sequence<uintE> unionEdge(uintE u,uintE v) {
    if(v_data.A[u].degree > v_data.A[v].degree) return unionEdge(v,u);
    auto fil = [&](uintE& t) {return existEdge(v,t);};
    auto entries = v_data.A[u].neighbors.entries();
    sequence<uintE> s = sequence<uintE>(v_data.A[u].degree,[&](uintE i){return std::get<0>(entries[i]);});
    return pbbslib::filter(s,fil);
  }

  void del() {
    deletion_fn();
  }

  vertex get_vertex(uintE i) {
    return vertex(v_data.A[i]);
  }

  void resizeVDATA(size_t s) {
    size_t nC = 1 << ((size_t)(ceil(log2(s))));
    dynamic_vertex_data* nA = pbbslib::new_array_no_init<dynamic_vertex_data>(nC);
    par_for(0,s,2000,[&](size_t i){nA[i] = v_data.A[i];});
    pbbslib::free_array(v_data.A);
    v_data.A = nA;
    v_data.capacity = nC;
    v_data.alloc = true;
  }

  void adjustVDATA(size_t amount) {
    size_t cur = v_data.capacity;
    if(amount != 0 && ceil(log2(amount)) == ceil(log2(cur))) {
      return;
    }
    resizeVDATA(amount);
  }

  void resizeNEIGHBORS(uintE u,size_t amount) {
    sparse_table<uintE,bool,hash_uintE> res = make_sparse_table<uintE,bool,hash_uintE>(amount,std::make_tuple(UINT_E_MAX,false),hash_uintE());
    auto entries = v_data.A[u].neighbors.entries();
    par_for(0,entries.size(),1,[&](size_t i) {
      uintE cur = std::get<0>(entries[i]);
      res.insert(std::make_tuple(cur,true));
    });
    v_data.A[u].neighbors = res;
  }

  void adjustNEIGHBORS(uintE u,size_t amount) {
    size_t cur = v_data.A[u].neighbors.m;
    if(amount != 0 && ceil(log2(amount)) == ceil(log2(cur))) {
      return;
    }
    resizeNEIGHBORS(u,amount);
  }

  void batchAddVertices(sequence<uintE> & vertices) {
    size_t size = vertices.size();

    pbbs::maxm<uintE> monoidm = pbbs::maxm<uintE>();
    uintE ma = pbbs::reduce(vertices,monoidm);

    adjustVDATA(std::max((size_t)ma,v_data.capacity));

    pbbs::sequence<uintE> existingVertices = pbbs::filter(vertices,[&](uintE i){return !existVertices.find(i,false);});
    size_t sumV = existingVertices.size();

    par_for(0,sumV,1,[&](size_t i) {
      uintE id = existingVertices[i];
      sparse_table<uintE,bool,hash_uintE> tmp = make_sparse_table<uintE,bool,hash_uintE>(INIT_VALUE_FOR_SIZE,std::make_tuple(UINT_E_MAX,false),hash_uintE());

      v_data.A[id].neighbors = tmp;
      v_data.A[id].degree = 0;
    });
    this -> n += sumV;
  }

  void addVertex(uintE id) {
    if(existVertices.find(id,false)) return;

    // TODO: Same issue with resizing as above.
    if(id >= v_data.size) v_data.resize((1 << (size_t)ceil(log2(id)) + 1) - v_data.size);

    sparse_table<uintE,bool,hash_uintE> tmp = make_sparse_table<uintE,bool,hash_uintE>(INIT_VALUE_FOR_SIZE,std::make_tuple(UINT_E_MAX,false),hash_uintE());
    
    v_data.A[id].neighbors = tmp;
    v_data.A[id].degree = 0;

    this -> n++;
    v_data.size++;
  }

  bool existEdge(uintE v,uintE u) {
    return v_data.A[v].neighbors.find(u,false);
  }

  void batchAddEdges(uintE u,sequence<uintE> edges) {
    pbbs::sequence<uintE> ds = pbbs::filter(edges,[&](uintE i){return !existEdge(u,i);});
    uintE sumW = ds.size();

    this -> m += sumW;
    adjustNEIGHBORS(u,v_data.A[u].degree + sumW);

    v_data.A[u].degree += ds.size();
    par_for(0,ds.size(),1,[&](size_t i) {
      uintE v = ds[i];
      adjustNEIGHBORS(v,v_data.A[u].degree + 1);

      v_data.A[v].neighbors.insert(std::make_tuple(u,true));
      v_data.A[u].neighbors.insert(std::make_tuple(v,true));

      v_data.A[v].degree++;
    });

  }

  void addEdge(uintE v,uintE u) {
    if(existEdge(u,v)) return;
    this -> m++;
    adjustNEIGHBORS(v,v_data.A[v].degree + 1);
    adjustNEIGHBORS(u,v_data.A[u].degree + 1);
    v_data.A[v].neighbors.insert(std::make_tuple(u,true));
    v_data.A[u].neighbors.insert(std::make_tuple(v,true));
  }

  void removeVertex(uintE id) {
    if(existVertices.find(id,false)) return;

    if(v_data.size - 1 < v_data.capacity >> 1) {
      v_data.resize((v_data.capacity >> 1) - v_data.size);
    }
    v_data.A[id].neighbors.del();
    v_data.A[id].degree = 0;
  }

  void batchRemoveVertices(sequence<uintE> & vertices) {
    size_t size = vertices.size();

    pbbs::sequence<uintE> existingVertices = pbbs::filter(vertices,[&](uintE i){return existVertices.find(i,false);});
    size_t sumV = existingVertices.size();
    par_for(0,sumV,1,[&](size_t i) {

      uintE id = existingVertices[i];

      v_data.A[id].neighbors.del();
      v_data.A[id].degree = 0;
    });

    if(v_data.size < v_data.capacity >> 1) {
      v_data.resize(v_data.capacity >> 1);
    }

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
      v_data.A[v].neighbors.erase(u);
      v_data.A[u].neighbors.erase(v);

      adjustNEIGHBORS(v,v_data.A[v].degree);
    });

    adjustNEIGHBORS(u,v_data.A[u].degree);


  }


  void removeEdge(uintE u,uintE v) {
    if(!existEdge(u,v)) return;
    this -> m--;

    v_data.A[u].degree--;
    v_data.A[v].degree--;
    v_data.A[v].neighbors.erase(u);
    v_data.A[u].neighbors.erase(v);
  }
};

template <template <class W> class vertex_type, class W>
dynamic_symmetric_graph<dynamic_symmetric_vertex,W> createEmptyDynamicSymmetricGraph() {

  std::function<void()> doNothing = []() {};
  pbbslib::dyn_arr<dynamic_vertex_data> _vData = pbbslib::dyn_arr<dynamic_vertex_data>(1);
  sparse_table<uintE,bool,hash_uintE> tmp =  make_sparse_table<uintE,bool,hash_uintE>(INIT_VALUE_FOR_SIZE,std::make_tuple(UINT_E_MAX,false),hash_uintE());

  return dynamic_symmetric_graph<dynamic_symmetric_vertex,W>(_vData,0,0,doNothing,tmp);
}


template <template <class W> class vertex_type, class W>
dynamic_symmetric_graph<dynamic_symmetric_vertex,W> dynamifyDSG(symmetric_graph<symmetric_vertex,pbbs::empty> G) {

  dynamic_symmetric_graph<dynamic_symmetric_vertex,W> dsg = createEmptyDynamicSymmetricGraph<dynamic_symmetric_vertex,W>();

  pbbs::sequence<uintE> v = pbbs::sequence<uintE>(G.n,[&](size_t i){return i;});
  dsg.batchAddVertices(v);

  for(int i = 0;i < G.n;i++) {
    pbbs::sequence<uintE> s = pbbs::sequence<uintE>(G.get_vertex(i).degree,[&](size_t j){return G.get_vertex(i).getOutNeighbor(j);});
    dsg.batchAddEdges(i,s);
  }
  
  return dsg;
}