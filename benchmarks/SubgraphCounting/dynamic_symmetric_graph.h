#pragma once

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
#include "ligra/pbbslib/sparse_table.h"
#include "ligra/pbbslib/dyn_arr.h"

// The initial sizes of dynamic arrays
size_t INIT_DYN_GRAPH_EDGE_SIZE = 4;

/**
 * Wrapper for hashing unsigned long long int
 */
struct hash_uintE {
  inline size_t operator () (const uintE & a) {return pbbs::hash64_2(a);}
};

/**
 * Struct of vertices. Stores the necessary information to recreate the node.
 */
struct dynamic_vertex_data {
  uintE degree;  // vertex degree
  sparse_table<uintE,bool,hash_uintE> neighbors; // Sparse_table
  uintE stored; // how many nodes are actually stored in the sparse table (some might be turned off)
};

/**
 * Wrapper of a vertex that has a lot of auxillary functions built in.
 * To invoke it, just pass in a dynamic_vertex_data.
 */
template <class W>
struct dynamic_symmetric_vertex {
  using vertex = dynamic_symmetric_vertex<W>;
  using edge_type = std::tuple<uintE, W>;
  using edges_type = sparse_table<uintE,bool,hash_uintE>;

  uintE degree;
  uintE stored;
  edges_type* neighbors = nullptr;

  dynamic_symmetric_vertex(dynamic_vertex_data& vdata) {
    neighbors = &vdata.neighbors;
    degree = vdata.degree;
    stored = vdata.stored;
  }

  // The function that should be called to delete the struct
  void clear() { neighbors -> del(); }
};

/**
* The struct that stores all the edges and vertices in the dynamic symmetric graph.
* It supports dynamic operations(adding and removing edges and vertices), and the graph is undirected(symmetric).
*
* To add vertices, use the batchAddVertices function, and to remove vertices, use the batchRemoveVertices function
* To add edges, use the batchAddEdges function, and to remove edges, use the batchRemoveEdges function
*
* @param W the template parameter that determins the type of the weight on the graph
* @param vertex_type the type of vertex data that should be used(it should usually just be dynamic_symmetric_vertex)
*/
template <template <class W> class vertex_type, class W>
struct dynamic_symmetric_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  using edge_type = typename vertex::edge_type;

  pbbslib::dyn_arr<dynamic_vertex_data> v_data; // List of the vertices, indexed by their vertex id.

  /* number of vertices in G */
  size_t n;
  /* number of edges in G */
  size_t m;

  // Boolean array that keeps track if a vertex exists
  pbbslib::dyn_arr<bool> existVertices;

  /* called to delete the graph */
  std::function<void()> deletion_fn;

  /**
  * The constructor for the dynamic symmetric graph.
  * This should be invoked to create the graph
  *
  * @param _v_data the list of existing vertices (can be empty if you are trying to create an empty graph)
  * @param _n the number of vertices in the graph
  * @param _m the number of edges in the graph
  * @param _deletion_fn the function that will be invoked once the struct is to be destructed
  * @param _existVertices this stores whether nodes are actually initialized in the v_data list.
  */
  dynamic_symmetric_graph(pbbslib::dyn_arr<dynamic_vertex_data> & _v_data, size_t _n, size_t _m,
    std::function<void()> _deletion_fn,pbbslib::dyn_arr<bool> & _existVertices)
      : n(_n),
        m(_m),
        deletion_fn(_deletion_fn) {
          v_data = pbbslib::dyn_arr<dynamic_vertex_data>(_v_data.size);
          v_data.size = _v_data.size;
          par_for(0,_v_data.size,1,[&](size_t i) {
            v_data.A[i] = _v_data.A[i];
          });
          existVertices = pbbslib::dyn_arr<bool>(_existVertices.size);
          existVertices.size = _existVertices.size;
          par_for(0,_existVertices.size,1,[&](size_t i) {
            existVertices.A[i] = _existVertices.A[i];
          });
  }

  /*
  * The function that should be when you want to delete the graph.
  */
  void del() {
    deletion_fn();
  }

  vertex get_vertex(uintE i) {
    return vertex(v_data.A[i]);
  }


  /**
  * Given the expected number of vertices, it would automatically adjust the v_data to an apporpriate capacity
  * 
  * @param amount The new expected number of vertices (or more accurately, the max ID expected to be stored in the graph)
  */
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
    par_for(0,nC,1,[&](size_t i){
      nAA[i] = i < existVertices.size ? existVertices.A[i] : false;
    });
    pbbslib::free_array(existVertices.A);
    existVertices.A = nAA;
    existVertices.capacity = nC;
  }



  /**
  * It would automatically adjust a node's data structures to an appropriate amount (it can both resize up and down)
  * 
  * @param u the index of the node that you want to be adjusted
  */
  void adjustNeighbors(uintE u) {
    size_t amount = std::max((size_t)(v_data.A[u].stored),INIT_DYN_GRAPH_EDGE_SIZE);
    size_t cur = v_data.A[u].neighbors.m;
    if((amount << 1) <= cur && (amount << 2) >= cur) return;
 
    auto entries = v_data.A[u].neighbors.entries();

    v_data.A[u].neighbors = make_sparse_table<uintE,bool,hash_uintE>(2 * amount,std::make_tuple(UINT_E_MAX,false),hash_uintE());
    par_for(0,entries.size(),1,[&](size_t i) {
      uintE cur = std::get<0>(entries[i]);
      if(std::get<1>(entries[i])) v_data.A[u].neighbors.insert(std::make_tuple(cur,true));
    });
    v_data.A[u].stored = v_data.A[u].degree;
  }

  /**
  * Adds vertices in parallel into the graph
  * 
  * @param vertices the sequence of vertices that you wish to be added into the graph
  */
  void batchAddVertices(sequence<uintE> & vertices) {

    size_t size = vertices.size();
    uintE ma = pbbs::reduce(vertices,pbbs::maxm<uintE>());
    adjustVdata(std::max((size_t)ma,v_data.capacity));


    pbbs::sequence<uintE> existingVertices = pbbs::filter(vertices,[&](uintE i){
      if (i < existVertices.size) return !existVertices.A[i];
      return true;
    });

    // cout << existVertices.capacity << " " << ma << endl;

    // cout << existingVertices.size() << endl;
    size_t sumV = existingVertices.size();

    par_for(0,sumV,1,[&](size_t i) {
      uintE id = existingVertices[i];
      existVertices.A[id] = true;

      v_data.A[id].neighbors = make_sparse_table<uintE,bool,hash_uintE>(INIT_DYN_GRAPH_EDGE_SIZE,std::make_tuple(UINT_E_MAX,false),hash_uintE());
      v_data.A[id].degree = 0;
      v_data.A[id].stored = 0;
    });
    n = std::max(this -> n,(size_t)ma + 1);
    v_data.size = n;
    existVertices.size = n;
  }


  // Check if an edge exists
  bool existEdge(uintE v,uintE u) {

    if (v < n && existVertices.A[v] && u < n && existVertices.A[u]) {
      return v_data.A[v].neighbors.find(u,false);
    }
    return false;
  }

  /**
  * Adds a series of edges to the graph in parallel that are centered around a specific node (faster than the other version)
  * The edges needs to be unique
  *
  * @param u the node that these edges are centered around
  * @param edges the list of nodes that will be u's neighbors
  */
  void batchAddEdges(uintE u,sequence<uintE> edges) {
    pbbs::sequence<uintE> ds = pbbs::filter(edges,[&](uintE i){return !existEdge(u,i);});
    uintE sumW = ds.size();

    this -> m += sumW;

    v_data.A[u].degree += ds.size();
    v_data.A[u].stored += ds.size();
    adjustNeighbors(u, sumW);

    par_for(0,ds.size(),1,[&](size_t i) {
      uintE v = ds[i];
      v_data.A[v].degree++;
      v_data.A[v].stored++;
      adjustNeighbors(v);

      v_data.A[v].neighbors.insert(std::make_tuple(u,true));
      v_data.A[u].neighbors.insert(std::make_tuple(v,true));
    });
  }


  /**
  * Adds a series of edges to the graph in parallel
  * The edges needs to be unique (i.e. (u,v) and (v,u) cannot be added at the same time)
  *
  * @param edges The sequence of edges that need to be removed
  */
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
    sequence<std::pair<uintE,uintE>> all = merge_sort(allS,[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a.first < b.first;});

    
    par_for(0,all.size(),[&](size_t i) {
      std::pair<uintE,uintE> cur = all[i];
      if(i != 0) {
        if(cur.first != all[i - 1].first) {
          v_data.A[all[i - 1].first].degree += i;
          v_data.A[all[i - 1].first].stored += i;
        }
      }
    });
    if(all.size() > 0) v_data.A[all[2 * ds.size() - 1].first].degree += 2 * ds.size();
    if(all.size() > 0) v_data.A[all[2 * ds.size() - 1].first].stored += 2 * ds.size();
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

    if(all.size() > 0) adjustNeighbors(all[0].first);

    par_for(0,all.size(),[&](size_t i) {
      std::pair<uintE,uintE> cur = all[i];
      v_data.A[cur.first].neighbors.insert(std::make_tuple(cur.second,true));
    });
  
  }

  /**
  * Removes vertices in parallel into the graph
  * 
  * @param vertices the sequence of vertices that you wish to be removed from the graph
  */
  void batchRemoveVertices(sequence<uintE> & vertices) {
    size_t size = vertices.size();

    pbbs::sequence<uintE> existingVertices = pbbs::filter(vertices,[&](uintE i){
      if (i < n) return existVertices.A[i];
      return false;
    });
    size_t sumV = existingVertices.size();
    par_for(0,sumV,1,[&](size_t i) {

      uintE id = existingVertices[i];
      existVertices.A[id] = false;

      v_data.A[id].neighbors.del();
      v_data.A[id].degree = 0;
    });

    adjustVdata(v_data.size);
  }

  /**
  * Removes a series of edges to the graph in parallel that are centered around a specific node (faster than the other version)
  * The edges needs to be unique
  *
  * @param u the node that these edges are centered around
  * @param edges the list of nodes that will cease to be u's neighbor
  */
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

  /**
  * Removes a series of edges to the graph in parallel
  * The edges needs to be unique (i.e. (u,v) and (v,u) cannot be added at the same time)
  *
  * @param edges The sequence of edges that need to be removed
  */
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
    sequence<std::pair<uintE,uintE>> all = merge_sort(allS,[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a.first < b.first; });
    
    par_for(0,all.size(),[&](size_t i) {
      std::pair<uintE,uintE> cur = all[i];
      v_data.A[cur.first].neighbors.change(cur.second,false);
      if(i != 0) {
        if(cur.first != all[i - 1].first) {
          v_data.A[all[i - 1].first].degree -= i;
        }
      }
    });
    if(all.size() > 0) v_data.A[all[2 * ds.size() - 1].first].degree -= 2 * ds.size();
    if(all.size() > 0) adjustNeighbors(all[0].first);

    par_for(0,all.size(),[&](size_t i) {
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

/**
* It returns an empty dynamic symmetric graph
* 
* @param W the template parameter that determins the type of the weight on the graph
* @param vertex_type the type of vertex data that should be used(it should usually just be dynamic_symmetric_vertex)
*/
template <template <class W> class vertex_type, class W>
dynamic_symmetric_graph<dynamic_symmetric_vertex,W> createEmptyDynamicSymmetricGraph() {

  std::function<void()> doNothing = []() {};
  pbbslib::dyn_arr<dynamic_vertex_data> _vData = pbbslib::dyn_arr<dynamic_vertex_data>(1);
  pbbslib::dyn_arr<bool> tmp = pbbslib::dyn_arr<bool>(1);
  tmp.A[0] = false;
  return dynamic_symmetric_graph<dynamic_symmetric_vertex,W>(_vData,0,0,doNothing,tmp);
}

/**
* It converts a standard gbbs graph to the dynamic symmetric version
* 
* @param W the template parameter that determins the type of the weight on the graph
* @param vertex_type the type of vertex data that should be used(it should usually just be dynamic_symmetric_vertex)
* @param G the original static gbbs Graph
*/
template <template <class W> class vertex_type, class W, class Graph>
dynamic_symmetric_graph<dynamic_symmetric_vertex,W> dynamifyDSG(Graph G) {

  dynamic_symmetric_graph<dynamic_symmetric_vertex,W> dsg = createEmptyDynamicSymmetricGraph<dynamic_symmetric_vertex,W>();

  pbbs::sequence<uintE> v = pbbs::sequence<uintE>(G.n,[&](size_t i){return i;});
  dsg.batchAddVertices(v);

  for(size_t i = 0; i < G.n;i++) {
    pbbs::sequence<uintE> s = pbbs::sequence<uintE>(G.get_vertex(i).degree,[&](size_t j){return G.get_vertex(i).getOutNeighbor(j);});
    dsg.batchAddEdges(i,s);
  }
  
  return dsg;
}