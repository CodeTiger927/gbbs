#include <string>

#include <stdlib.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <cmath>

#include "pbbslib/sequence_ops.h"

#include "ligra/pbbslib/sparse_table.h"

#include "ligra/pbbslib/dyn_arr.h"



#define INIT_VALUE_FOR_SIZE 100



// General Comment and TODO (From Alex):
// After coding the Symmetric Graph Class. I cannot help but notice that the vertex class is
// pretty useless since all of its functionality can be easily replaced by directly using
// vertex data, and due to the nature of sparse_table, I have to create a entries dyn_arr
// every time we initialize a vertex, which is like O(N) at first, so maybe I should somehow
// remove vertex or integrate the entries dyn_arr into vertex_data.



template <template <typename W> class vertex, class W, class F>
inline void mapNghs(vertex<W>* v, uintE vtx_id, std::tuple<uintE, W>* nghs,
                    uintE d, F& f, bool parallel) {
  par_for(0, d, pbbslib::kSequentialForThreshold, [&] (size_t j) {
    uintE ngh = v->getOutNeighbor(j);
    f(vtx_id, ngh, v->getOutWeight(j));
  }, parallel);
}





std::function<void()> get_deletion_fn(void*, void*);
std::function<void()> get_deletion_fn(void*, void*, void*);
std::function<void()> get_deletion_fn(void*, void*, void*, void*);

/* Compressed Sparse Row (CSR) based representation for symmetric graphs.
 * Takes two template parameters:
 * 1) vertex_type: vertex template, parametrized by the weight type associated with each edge
 * 2) W: the weight template
 * The graph is represented as an array of edges of type vertex_type::edge_type,
 * which is just a tuple<uintE, W>.*/


// currently this is using CSR so I need to change this into sparse table




// We need to remove this in one of the files to avoid redefinition (h-index2.h)
struct hash_uintE {
  inline size_t operator () (const uintE & a) {return pbbs::hash64_2(a);}
};



struct dynamic_vertex_data {
  uintE degree;  // vertex degree
  sparse_table<uintE,bool,hash_uintE> neighbors;
  pbbslib::dyn_arr<std::tuple<uintE, uintE>> entries;
};


// This is just the vertex class itself


template <class W>
struct dynamic_symmetric_vertex {
  using vertex = dynamic_symmetric_vertex<W>;
  using edge_type = std::tuple<uintE, W>;
  using edges_type = sparse_table<uintE,bool,hash_uintE>;

  uintE degree;
  edges_type neighbors = make_sparse_table<uintE,bool,hash_uintE>(INIT_VALUE_FOR_SIZE,std::make_tuple(UINT_E_MAX,false),hash_uintE());
  pbbslib::dyn_arr<edge_type> entries;

  dynamic_symmetric_vertex(dynamic_vertex_data& vdata) {
    neighbors = vdata.neighbors;
    degree = vdata.degree;
    entries = vdata.entries;

    // sequence<std::tuple<uintE,bool>> s = neighbors.entries();
    // par_for(0,s.size(),1,[&](size_t i) {
    //   if(std::get<1>(s[i])) this -> entries.push_back(std::make_tuple(std::get<0>(s[i]),true));
    // });
  }


  // Testing, can be deleted later
  dynamic_symmetric_vertex() {
    degree = 0;
  }


  void clear() { neighbors.del(); }


  edges_type getInNeighbors() { return neighbors; }
  edges_type getOutNeighbors() { return neighbors; }
  uintE getInNeighbor(uintE j) { return std::get<0>(this -> entries.A[j]); }
  uintE getOutNeighbor(uintE j) { return std::get<0>(this -> entries.A[j]); }
  W getInWeight(uintE j) { return std::get<1>(this -> entries.A[j]); }
  W getOutWeight(uintE j) { return std::get<1>(this -> entries.A[j]); }

  uintE getInDegree() { return degree; }
  uintE getOutDegree() { return degree; }
  uintE getInVirtualDegree() { return degree; }
  uintE getOutVirtualDegree() { return degree; }

  template <class F>
  inline void mapOutNgh(uintE vtx_id, F& f, bool parallel = true) {
    mapNghs<dynamic_symmetric_vertex, W, F>(this, vtx_id, getOutNeighbors(),
                                               getOutDegree(), f, parallel);
  }

  template <class F>
  inline void mapInNgh(uintE vtx_id, F& f, bool parallel = true) {
    mapNghs<dynamic_symmetric_vertex, W, F>(this, vtx_id, getInNeighbors(),
                                               getOutDegree(), f, parallel);
  }

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

  /* called to delete the graph */
  std::function<void()> deletion_fn;

  dynamic_symmetric_graph(pbbslib::dyn_arr<dynamic_vertex_data> v_data, size_t n, size_t m,
    std::function<void()> _deletion_fn)
      : v_data(v_data),
        n(n),
        m(m),
        deletion_fn(_deletion_fn) {
  }


  sequence<uintE> unionEdge(uintE u,uintE v) {
    if(v_data.A[u].entries.size > v_data.A[v].entries.size) return unionEdge(v,u);
    auto fil = [&](uintE& t) {return existEdge(v,t);};
    sequence<uintE> s = sequence<uintE>(v_data.A[u].entries.size,[&](uintE i){return std::get<0>(v_data.A[u].entries.A[i]);});
    return pbbslib::filter(s,fil);
  }


  void del() {
    deletion_fn();
  }

  vertex get_vertex(uintE i) {
    return vertex(v_data.A[i]);
  }

// #ifndef TWOSOCKETNVM
//   vertex get_vertex(uintE i) {
//     return vertex(v_data[i]);
//   }
// #else
//   vertex get_vertex(uintE i) {
//     if (numanode() == 0) {
//       return vertex(v_data[i]);
//     } else {
//       return vertex(v_data[i]);
//     }
//   }
// #endif

  /* degree must be <= old_degree */
  void decrease_vertex_degree(uintE id, uintE degree) {
    assert(degree <= v_data.A[id].degree);
    v_data.A[id].degree = degree;
  }

  void zero_vertex_degree(uintE id) {
    decrease_vertex_degree(id, 0);
  }

  void batchAddVertexes(sequence<uintE> & vertexes) {
    size_t size = vertexes.size();

    uintE m;
    // From Alex: Hmm so basically I dont think this is
    // supposed to be working since what if multiple parallel
    // processors try to change m at the same time?
    // But it seems to work fine when I test it? Idk
    m = 0;
    par_for(0,size,1,[&](size_t i) {
      if(vertexes[i] > m) m = vertexes[i];
    });

    if(m >= v_data.size) {
      v_data.resize(1 << (int)floor(log(m)) + 1);
    }    
    par_for(0,size,1,[&](size_t i) {
      uintE id = vertexes[i];
      sparse_table<uintE,bool,hash_uintE> tmp = make_sparse_table<uintE,bool,hash_uintE>(INIT_VALUE_FOR_SIZE,std::make_tuple(UINT_E_MAX,false),hash_uintE());
      pbbslib::dyn_arr<std::tuple<uintE, uintE>> entries = pbbslib::dyn_arr<std::tuple<uintE, uintE>>(0);
    
      v_data.A[id].neighbors = tmp;
      v_data.A[id].degree = 0;
      v_data.A[id].entries = entries;
    });

    this -> n += size;
  }

  // TODO: Put in a wrapper function that does batch-parallel vertex and 
  // edge additions/deletions. It should take as input an array or function, and a size,
  // and it will add and delete all vertices/edges in that input in parallel.
  void addVertex(uintE id) {

    if(id >= v_data.size) v_data.resize(1 << (int)floor(log(id)) + 1);


    // TODO: Do not use a constant like 1000 to initialize. You should initialize
    // to a reasonable value that is either a constant variable or an input, 
    // and you'll need a way to resize.
    sparse_table<uintE,bool,hash_uintE> tmp = make_sparse_table<uintE,bool,hash_uintE>(INIT_VALUE_FOR_SIZE,std::make_tuple(UINT_E_MAX,false),hash_uintE());
    pbbslib::dyn_arr<std::tuple<uintE, uintE>> entries = pbbslib::dyn_arr<std::tuple<uintE, uintE>>(0);
    
    v_data.A[id].neighbors = tmp;
    v_data.A[id].degree = 0;
    v_data.A[id].entries = entries;

    this -> n++;
  }

  bool existEdge(uintE v,uintE u) {
    return v_data.A[v].neighbors.find(u,false);
  }

  void _checkSize(uintE v,uintE many) {
    if(v_data.A[v].entries.size > v_data.A[v].entries.capacity - 10 - many) {
      v_data.A[v].entries.resize(1 << (int)floor(v_data.A[v].entries.size + 10 + many) + 1);
    }
    if(v_data.A[v].degree > v_data.A[v].neighbors.m - 10 - many) {
      sparse_table<uintE,bool,hash_uintE> tmp = make_sparse_table<uintE,bool,hash_uintE>(1 << (int)floor(v_data.A[v].neighbors.m + many + 10) + 1,std::make_tuple(UINT_E_MAX,false),hash_uintE());
      par_for(0,v_data.A[v].entries.size,1,[&](size_t i) {
        tmp.insert(std::make_tuple(std::get<0>(v_data.A[v].entries.A[i]),true));
      });
      v_data.A[v].neighbors = tmp;
    }
  }

  void batchAddEdges(uintE u,sequence<uintE> edges) {

    this -> m += edges.size();
    _checkSize(u,edges.size());
    v_data.A[u].degree += edges.size();
    par_for(0,edges.size(),1,[&](size_t i) {
      uintE v = edges[i];
      _checkSize(v,1);

      v_data.A[v].neighbors.insert(std::make_tuple(u,true));
      v_data.A[v].entries.push_back(std::make_tuple(u,0));
      v_data.A[u].entries.push_back(std::make_tuple(v,0));
      v_data.A[u].neighbors.insert(std::make_tuple(v,true));

      v_data.A[v].degree++;
    });
  }

  void addEdge(uintE v,uintE u) {
    // TODO: Add something to prevent non-existant vertexes to connect each other. Or if an already existing edge
    this -> m++;
    _checkSize(v,1);
    _checkSize(u,1);
    v_data.A[v].neighbors.insert(std::make_tuple(u,true));
    v_data.A[v].entries.push_back(std::make_tuple(u,0));
    v_data.A[u].neighbors.insert(std::make_tuple(v,true));
    v_data.A[u].entries.push_back(std::make_tuple(v,0));
  }


  template <class F>
  void _edges(F f, bool parallel_inner_map = true) {
    parallel_for(0, n, [&](size_t i) {
      get_vertex(i).mapOutNgh(i, f, parallel_inner_map);
    }, 1);
  }

  // F : edge -> edge
  template <class F>
  void alter_edges(F f, bool parallel_inner_map = true) {
    abort(); /* unimplemented for CSR */
  }





  pbbs::sequence<std::tuple<uintE, uintE, W>> edges() {
    using g_edge = std::tuple<uintE, uintE, W>;
    auto degs = pbbs::sequence<size_t>(n, [&] (size_t i) {
      return get_vertex(i).getOutDegree();
    });
    size_t sum_degs = pbbslib::scan_add_inplace(degs.slice());
    assert(sum_degs == m);
    auto edges = pbbs::sequence<g_edge>(sum_degs);
    parallel_for(0, n, [&](size_t i) {
      size_t k = degs[i];
      auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
        edges[k++] = std::make_tuple(u, v, wgh);
      };
      get_vertex(i).mapOutNgh(i, map_f, false);
    }, 1);
    return edges;
  }


};



// // Mutates (sorts) the underlying array A containing a black-box description of
// // an edge of typename A::value_type. The caller provides functions GetU, GetV, and GetW
// // which extract the u, v, and weight of a (u,v,w) edge respective (if the edge
// // is a std::tuple<uinte, uintE, W> this is just get<0>, ..<1>, ..<2>
// // respectively.
// // e.g.:
// //   using edge = std::tuple<uintE, uintE, W>;
// //   auto get_u = [&] (const edge& e) { return std::get<0>(e); };
// //   auto get_v = [&] (const edge& e) { return std::get<1>(e); };
// //   auto get_w = [&] (const edge& e) { return std::get<2>(e); };
// //   auto G = sym_graph_from_edges<W>(coo1, get_u, get_v, get_w, 10, false);
// template <class W, class EdgeSeq, class GetU, class GetV, class GetW>
// inline dynamic_symmetric_graph<dynamic_symmetric_vertex, W> sym_graph_from_edges(
//     EdgeSeq& A,
//     size_t n,
//     GetU&& get_u,
//     GetV&& get_v,
//     GetW&& get_w,
//     bool is_sorted = false) {
//   using vertex = dynamic_symmetric_vertex<W>;
//   using edge_type = typename vertex::edge_type;
//   size_t m = A.size();

//   if (m == 0) {
//     if (n == 0) {
//       std::function<void()> del = []() {};
//       return dynamic_symmetric_graph<dynamic_symmetric_vertex, W>(nullptr, 0, 0, del, nullptr);
//     } else {
//       auto v_data = pbbs::new_array_no_init<dynamic_vertex_data>(n);
//       par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
//         v_data[i].offset = 0;
//         v_data[i].degree = 0;
//       });
//       std::function<void()> del = get_deletion_fn(v_data, nullptr);
//       return dynamic_symmetric_graph<dynamic_symmetric_vertex, W>(v_data, n, 0, del, nullptr);
//     }
//   }

//   if (!is_sorted) {
//     size_t bits = pbbslib::log2_up(n);
//     pbbslib::integer_sort_inplace(A.slice(), get_u, bits);
//   }

//   auto starts = sequence<uintT>(n+1, (uintT) 0);

//   using neighbor = std::tuple<uintE, W>;
//   auto edges = sequence<neighbor>(m, [&](size_t i) {
//     // Fuse loops over edges (check if this helps)
//     if (i == 0 || (get_u(A[i]) != get_u(A[i - 1]))) {
//       starts[get_u(A[i])] = i;
//     }
//     if (i != (m-1)) {
//       uintE our_vtx = get_u(A[i]);
//       uintE next_vtx = get_u(A[i+1]);
//       if (our_vtx != next_vtx && (our_vtx + 1 != next_vtx)) {
//         par_for(our_vtx+1, next_vtx, pbbslib::kSequentialForThreshold, [&] (size_t k) {
//           starts[k] = i+1;
//         });
//       }
//     }
//     if (i == (m-1)) { /* last edge */
//       starts[get_u(A[i]) + 1] = m;
//     }
//     return std::make_tuple(get_v(A[i]), get_w(A[i]));
//   });

// //  auto copy_f = [&] (const uintT& u, const uintT& v) -> uintT {
// //    if (v == 0) { return u; }
// //    return v;
// //  };
// //  auto copy_m = pbbs::make_monoid(copy_f, (uintT)0);
// //  pbbs::scan_inplace(starts.slice(), copy_m, pbbs::fl_inplace);
// //  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
// //    uintT o = starts[i];
// //    size_t degree = ((i == n - 1) ? m : starts[i + 1]) - o;
// //    v[i].degree = degree;
// //    v[i].neighbors = ((std::tuple<uintE, W>*)(edges.begin() + o));
// //    /* sort each neighbor list if needed */
// //  });
// //  return graph<V>(v, n, m, get_deletion_fn(v, edges.to_array()));

//   auto v_data = pbbs::new_array_no_init<dynamic_vertex_data>(n);
//   par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
//     uintT o = starts[i];
//     v_data[i].offset = o;
//     v_data[i].degree = (uintE)(((i == (n-1)) ? m : starts[i+1]) - o);
//   });
//   auto new_edge_arr = edges.to_array();
//   return dynamic_symmetric_graph<dynamic_symmetric_vertex, W>(v_data, n, m, get_deletion_fn(v_data, new_edge_arr), (edge_type*)new_edge_arr);
// }

// template <class W>
// inline dynamic_symmetric_graph<dynamic_symmetric_vertex, W> sym_graph_from_edges(
//     pbbs::sequence<std::tuple<uintE, uintE, W>>& A,
//     size_t n,
//     bool is_sorted = false) {
//   using edge = std::tuple<uintE, uintE, W>;
//   auto get_u = [&] (const edge& e) { return std::get<0>(e); };
//   auto get_v = [&] (const edge& e) { return std::get<1>(e); };
//   auto get_w = [&] (const edge& e) { return std::get<2>(e); };
//   return sym_graph_from_edges<W>(A, n, get_u, get_v, get_w, is_sorted);
// }



template <template <class W> class vertex_type, class W>
dynamic_symmetric_graph<dynamic_symmetric_vertex,W> createEmptyDynamicSymmetricGraph() {

  std::function<void()> doNothing = []() {};
  pbbslib::dyn_arr<dynamic_vertex_data> _vData = pbbslib::dyn_arr<dynamic_vertex_data>(0);

  return dynamic_symmetric_graph<dynamic_symmetric_vertex,W>(_vData,0,0,doNothing);
}