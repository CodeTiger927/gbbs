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



// TODO: It's not really good practice to use a preprocessor directive to set
// a constant. If anything, you should use a global constant. In reality,
// this value should not be a constant at all -- it should be chosen dynamically
// based off of the graph G and the initial max degree or degree distribution
// of the graph.
uintE INIT_VALUE_FOR_SIZE = 10;



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


// We need to remove this in one of the files to avoid redefinition (h-index2.h)
struct hash_uintE {
  inline size_t operator () (const uintE & a) {return pbbs::hash64_2(a);}
};



struct dynamic_vertex_data {
  uintE degree;  // vertex degree
  sparse_table<uintE,uintE,hash_uintE> neighbors;
  pbbslib::dyn_arr<std::tuple<uintE,uintE>> entries;
};


// This is just the vertex class itself


template <class W>
struct dynamic_symmetric_vertex {
  using vertex = dynamic_symmetric_vertex<W>;
  using edge_type = std::tuple<uintE, W>;
  using edges_type = sparse_table<uintE,uintE,hash_uintE>;

  uintE degree;
  edges_type neighbors = make_sparse_table<uintE,uintE,hash_uintE>(INIT_VALUE_FOR_SIZE,std::make_tuple(UINT_E_MAX,UINT_E_MAX),hash_uintE());
  pbbslib::dyn_arr<edge_type> entries;

  dynamic_symmetric_vertex(dynamic_vertex_data& vdata) {
    neighbors = vdata.neighbors;
    degree = vdata.degree;
    entries = vdata.entries;
  }


  // Testing, can be deleted later
  dynamic_symmetric_vertex() {
    degree = 0;
  }


  void clear() { neighbors.del(); }


  edges_type getInNeighbors() { return neighbors; }
  edges_type getOutNeighbors() { return neighbors; }
  // TODO: I don't think you need to use "this" since it's unambiguous what entries is
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

  // TODO: Don't overload names; stylistically, we would use _v_data as the input
  // name, and set v_data(_v_data) or something along those lines.
  dynamic_symmetric_graph(pbbslib::dyn_arr<dynamic_vertex_data> &v_data, size_t n, size_t m,
    std::function<void()> _deletion_fn)
      : n(n),
        m(m),
        deletion_fn(_deletion_fn) {
          this -> v_data = pbbslib::dyn_arr<dynamic_vertex_data>(v_data.size);
          par_for(0,v_data.size,1,[&](size_t i) {
            this -> v_data.A[i] = v_data.A[i];
          });
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


  void batchAddVertices(sequence<uintE> & vertices) {
    size_t size = vertices.size();


    // Find max ew vertices that must be added
    pbbs::maxm<uintE> monoidm = pbbs::maxm<uintE>();
    uintE ma = pbbs::reduce(vertices,monoidm);

    if(ma >= v_data.capacity) {
      v_data.resize((1 << (size_t)ceil(log2(ma)) + 1) - v_data.size);
    }

    // TODO: Maybe it would be better here to use a delayed sequence, so we're not actually creating all of this space.
    // Alternatively, don't just use addm in the reduce -- iterate through v_data.A or something and write your
    // own monoid to retrieve whether it's allocated or not, and pass the monoid directly into reduce.
    pbbs::sequence<uintE> ds = pbbs::sequence<uintE>(size,[&](size_t i){return !v_data.A[vertices[i]].entries.alloc;});
    uintE sumV = pbbs::reduce(ds, pbbs::addm<uintE>());
    // TODO: What are you doing with v_data.size here? /Why are you overloading v_data.size to keep track of which
    // vertex_data entries have been allocated? Shouldn't the size be given when you do a batch addition to 
    // dyn_arr?
    v_data.size += sumV;

    par_for(0,size,1,[&](size_t i) {
      if(!v_data.A[vertices[i]].entries.alloc) {

        uintE id = vertices[i];
        sparse_table<uintE,uintE,hash_uintE> tmp = make_sparse_table<uintE,uintE,hash_uintE>(INIT_VALUE_FOR_SIZE,std::make_tuple(UINT_E_MAX,UINT_E_MAX),hash_uintE());
        pbbslib::dyn_arr<std::tuple<uintE, uintE>> entries = pbbslib::dyn_arr<std::tuple<uintE, uintE>>(1);
      
        v_data.A[id].neighbors = tmp;
        v_data.A[id].degree = 0;
        v_data.A[id].entries = entries;

      }
    });


    this -> n += sumV;

  }

  void addVertex(uintE id) {

    if(v_data.size > id && v_data.A[id].entries.alloc) return;

    // TODO: Same issue with resizing as above.
    if(id >= v_data.size) v_data.resize((1 << (size_t)ceil(log2(id)) + 1) - v_data.size);

    sparse_table<uintE,bool,hash_uintE> tmp = make_sparse_table<uintE,bool,hash_uintE>(INIT_VALUE_FOR_SIZE,std::make_tuple(UINT_E_MAX,false),hash_uintE());
    pbbslib::dyn_arr<std::tuple<uintE, uintE>> entries = pbbslib::dyn_arr<std::tuple<uintE, uintE>>(0);
    
    v_data.A[id].neighbors = tmp;
    v_data.A[id].degree = 0;
    v_data.A[id].entries = entries;

    this -> n++;
    v_data.size++;
  }

  bool existEdge(uintE v,uintE u) {
    return v_data.A[v].neighbors.find(u,UINT_E_MAX) != UINT_E_MAX;
  }

  // TODO: What is this meant to do? Why the -10?
  void _checkSize(uintE v,uintE many) {
    if(v_data.A[v].entries.size + many > v_data.A[v].entries.capacity) {
      v_data.A[v].entries.resize(many);
      if(v == 1) {
        cout << v_data.A[v].entries.capacity << endl;
      }
    }
    if(v_data.A[v].entries.size + many > v_data.A[v].neighbors.m) {
      sparse_table<uintE,uintE,hash_uintE> tmp = make_sparse_table<uintE,uintE,hash_uintE>(v_data.A[v].entries.size + many,std::make_tuple(UINT_E_MAX,UINT_E_MAX),hash_uintE());
      par_for(0,v_data.A[v].entries.size,1,[&](size_t i) {
        tmp.insert(std::make_tuple(std::get<0>(v_data.A[v].entries.A[i]),i));
      });
      v_data.A[v].neighbors = tmp;
    }
  }


  // CheckSize For deletion
  void __checkSize(uintE v,uintE many) {
    // TODO: Need to write an effective shrinking function to conserve memory
  }

  void batchAddEdges(uintE u,sequence<uintE> edges) {
    pbbs::sequence<uintE> ds = pbbs::filter(edges,[&](uintE i){return !existEdge(u,i);});
    uintE sumW = ds.size();

    this -> m += sumW;
    _checkSize(u,sumW);


    size_t ori = v_data.A[u].entries.size;
    v_data.A[u].entries.size += ds.size();

    v_data.A[u].degree += ds.size();
    par_for(0,ds.size(),1,[&](size_t i) {
      uintE v = ds[i];
      _checkSize(v,1);

        //cout << "No segmentation fault yet " << i << endl;
      v_data.A[v].neighbors.insert(std::make_tuple(u,v_data.A[v].entries.size));
      v_data.A[v].entries.push_back(std::make_tuple(u,0));
      v_data.A[u].entries.A[ori + i] = std::make_tuple(v,0);
      v_data.A[u].neighbors.insert(std::make_tuple(v,ori + i));

      v_data.A[v].degree++;
    });
  }



  // Btw this cannot be used in parallel.
  void addEdge(uintE v,uintE u) {
    // TODO: Add something to prevent non-existant vertices to connect each other. Or if an already existing edge
    //cout << "YES  " << u << " " << v << " " << existEdge(v,u) << endl;
    if(existEdge(u,v)) return;
    this -> m++;
    _checkSize(v,1);
    _checkSize(u,1);
    v_data.A[v].neighbors.insert(std::make_tuple(u,v_data.A[v].entries.size));
    v_data.A[v].entries.push_back(std::make_tuple(u,0));
    v_data.A[u].neighbors.insert(std::make_tuple(v,v_data.A[u].entries.size));
    v_data.A[u].entries.push_back(std::make_tuple(v,0));
  }

  // TODO, check if the vertex has degree, if so then cannot delete
  // Btw this cannot be used in parallel too.
  void removeVertex(uintE id) {
    if(v_data.size > id && !(v_data.A[id].entries.alloc)) return;

    if(v_data.size - 1 < v_data.capacity >> 1) {
      v_data.resize((v_data.capacity >> 1) - v_data.size);
    }
    v_data.A[id].neighbors.del();
    v_data.A[id].degree = 0;
    v_data.A[id].entries.del();
  }

  void batchRemoveVertices(sequence<uintE> & vertices) {
    size_t size = vertices.size();

    pbbs::sequence<uintE> ds = pbbs::sequence<uintE>(size,[&](size_t i){return v_data.A[vertices[i]].entries.alloc;});
    uintE sumV = pbbs::reduce(ds, pbbs::addm<uintE>());


    v_data.size -= sumV;

    par_for(0,size,1,[&](size_t i) {
      if(v_data.A[vertices[i]].entries.alloc) {

        uintE id = vertices[i];

        v_data.A[id].neighbors.del();
        v_data.A[id].degree = 0;
        v_data.A[id].entries.del();
      }
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

    size_t ori = v_data.A[u].entries.size;
    v_data.A[u].entries.size -= ds.size();
    v_data.A[u].degree -= ds.size();
    sparse_table<uintE,bool,hash_uintE> checkIfThere = make_sparse_table<uintE,bool,hash_uintE>(ds.size(),std::make_tuple(UINT_E_MAX,false),hash_uintE());

    par_for(0,ds.size(),1,[&](size_t i) {
      uintE v = ds[i];
      checkIfThere.insert(std::make_tuple(v,true));

      // Remove u from v
      uintE whereU = v_data.A[v].neighbors.find(u,UINT_E_MAX);
      v_data.A[v].entries.size--;
      v_data.A[v].degree--;
      v_data.A[v].entries.A[whereU] = v_data.A[v].entries.A[v_data.A[v].entries.size];
      v_data.A[v].neighbors.erase(u);

      if(whereU != v_data.A[v].entries.size) v_data.A[v].neighbors.erase(std::get<0>(v_data.A[v].entries.A[whereU]));

      if(whereU != v_data.A[v].entries.size) v_data.A[v].neighbors.insert(std::make_tuple(std::get<0>(v_data.A[v].entries.A[whereU]),whereU));

      __checkSize(v,1);

    });

    // A separate loop to avoid collisions and segmentation fault
    pbbs::sequence<uintE> notInLast = pbbs::filter(ds,[&](uintE i) {
      return (v_data.A[u].neighbors.find(i,UINT_E_MAX) < (ori - sumW));
    });
    pbbs::sequence<uintE> allLast = pbbs::sequence<uintE>(sumW,[&](size_t i){return (ori - i - 1);});
    pbbs::sequence<uintE> notInDS = pbbs::filter(allLast,[&](uintE i){return !checkIfThere.find(std::get<0>(v_data.A[u].entries.A[i]),false);});

    for(int i = 0;i < notInLast.size();i++) {
      // swap notInLast[i] and notInDS[i]
      uintE whereU = v_data.A[u].neighbors.find(notInLast[i],UINT_E_MAX);
      uintE whereV = notInDS[i];
      v_data.A[u].entries.A[whereU] = v_data.A[u].entries.A[whereV];
      v_data.A[u].neighbors.change(notInLast[i],UINT_E_MAX);
      v_data.A[u].neighbors.change(std::get<0>(v_data.A[u].entries.A[whereV]),whereU);
    }
    __checkSize(u,sumW);
  }

  void removeEdge(uintE u,uintE v) {
    if(!existEdge(u,v)) return;
    this -> m--;
    __checkSize(u,1);
    __checkSize(v,1);

    uintE whereV = v_data.A[u].neighbors.find(v,UINT_E_MAX);
    uintE whereU = v_data.A[v].neighbors.find(u,UINT_E_MAX);

    v_data.A[u].entries.size--;
    v_data.A[v].entries.size--;

    v_data.A[u].degree--;
    v_data.A[v].degree--;

    v_data.A[u].entries.A[whereV] = v_data.A[u].entries.A[v_data.A[u].entries.size];
    v_data.A[v].entries.A[whereU] = v_data.A[v].entries.A[v_data.A[v].entries.size];

    v_data.A[u].neighbors.erase(v);
    v_data.A[v].neighbors.erase(u);


    if(whereV != v_data.A[u].entries.size) v_data.A[u].neighbors.erase(std::get<0>(v_data.A[u].entries.A[whereV]));
    if(whereU != v_data.A[v].entries.size) v_data.A[v].neighbors.erase(std::get<0>(v_data.A[v].entries.A[whereU]));

    if(whereV != v_data.A[u].entries.size) v_data.A[u].neighbors.insert(std::make_tuple(std::get<0>(v_data.A[u].entries.A[whereV]),whereV));
    if(whereU != v_data.A[v].entries.size) v_data.A[v].neighbors.insert(std::make_tuple(std::get<0>(v_data.A[v].entries.A[whereU]),whereU));
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


template <template <class W> class vertex_type, class W>
dynamic_symmetric_graph<dynamic_symmetric_vertex,W> createEmptyDynamicSymmetricGraph() {

  std::function<void()> doNothing = []() {};
  pbbslib::dyn_arr<dynamic_vertex_data> _vData = pbbslib::dyn_arr<dynamic_vertex_data>(1);

  return dynamic_symmetric_graph<dynamic_symmetric_vertex,W>(_vData,0,0,doNothing);
}


template <template <class W> class vertex_type, class W>
dynamic_symmetric_graph<dynamic_symmetric_vertex,W> dynamifyDSG(symmetric_graph<symmetric_vertex,pbbs::empty> G) {
  std::function<void()> doNothing = []() {};
  pbbslib::dyn_arr<dynamic_vertex_data> _vData = pbbslib::dyn_arr<dynamic_vertex_data>(1);

  dynamic_symmetric_graph<dynamic_symmetric_vertex,W> dsg = dynamic_symmetric_graph<dynamic_symmetric_vertex,W>(_vData,0,0,doNothing);

  pbbs::sequence<uintE> v = pbbs::sequence<uintE>(G.n,[&](size_t i){return i;});
  dsg.batchAddVertices(v);

  for(int i = 0;i < G.n;i++) {
    pbbs::sequence<uintE> s = pbbs::sequence<uintE>(G.get_vertex(i).degree,[&](size_t j){return G.get_vertex(i).getOutNeighbor(j);});
    dsg.batchAddEdges(i,s);
  }
  
  return dsg;
}