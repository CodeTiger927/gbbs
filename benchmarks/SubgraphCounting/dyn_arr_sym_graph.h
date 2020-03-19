#include "ligra/pbbslib/dyn_arr.h"


//No need for vertex_data because it would store dyn_arr of neighbors anyways
struct dyn_arr_sym_vertex {
  uintE id;
  pbbslib::dyn_arr<uintE> neighbors; 
  uintE deg;
  //Doesn't care about edges because it has neighbors and weight doesn't matter

  //n - initial capacity of neighbors array
  //Will be resized if necessary
  dyn_arr_sym_vertex(uintE _id, size_t n):
    id(_id),
    deg(0),
    neighbors(pbbslib::dyn_arr<uintE>(n))
  {}

  void add_neighbor(uintE v) {
    deg++;
    neighbors.add(v);
  }

  void remove_neighbor(uintE v) {
    deg--;
    neighbors.erase(v);
  }

  //Do we need all this
  /*
  void clear() { pbbslib::free_array(neighbors); }

  edge_type* getInNeighbors() { return neighbors; }
  edge_type* getOutNeighbors() { return neighbors; }
  uintE getInNeighbor(uintE j) { return std::get<0>(neighbors[j]); }
  uintE getOutNeighbor(uintE j) { return std::get<0>(neighbors[j]); }
  W getInWeight(uintE j) { return std::get<1>(neighbors[j]); }
  W getOutWeight(uintE j) { return std::get<1>(neighbors[j]); }

  uintE getInDegree() { return degree; }
  uintE getOutDegree() { return degree; }
  uintE getInVirtualDegree() { return degree; }
  uintE getOutVirtualDegree() { return degree; }

  constexpr static uintE getInternalBlockSize() {
    return vertex_ops::kBlockSize;
  }
  uintE getNumInBlocks() {
    return pbbs::num_blocks(degree, vertex_ops::kBlockSize);
  }
  uintE getNumOutBlocks() { return getNumInBlocks(); }
  inline uintE in_block_degree(uintE block_num) {
    uintE block_start = block_num * vertex_ops::kBlockSize;
    uintE block_end = std::min(block_start + vertex_ops::kBlockSize, degree);
    return block_end - block_start;  // TODO: check
  }
  inline uintE out_block_degree(uintE block_num) {
    return in_block_degree(block_num);
  }

  void flipEdges() {}

  auto getInIter(uintE id) -> vertex_ops::iter<W> {
    return vertex_ops::get_iter(getInNeighbors(), getInDegree());
  }
  auto getOutIter(uintE id) -> vertex_ops::iter<W> {
    return vertex_ops::get_iter(getOutNeighbors(), getOutDegree());
  }

  inline size_t intersect(symmetric_vertex<W>* other, long our_id,
                          long other_id) {
    return intersection::intersect(this, other, our_id, other_id);
  }

  template <class F>
  inline size_t intersect_f(symmetric_vertex<W>* other, long our_id,
                            long other_id, const F& f) {
    return intersection::intersect_f(this, other, our_id, other_id, f);
  }

  template <class F>
  inline size_t intersect_f_par(symmetric_vertex<W>* other, long our_id,
                            long other_id, const F& f) {
    return intersection::intersect_f_par(this, other, our_id, other_id, f);
  }


  template <class VS, class F, class G>
  inline void decodeInNghBreakEarly(uintE vtx_id, VS& vertexSubset, F& f, G& g,
                                    bool parallel = 0) {
    vertex_ops::decodeNghsBreakEarly<symmetric_vertex, W, F, G, VS>(
        this, vtx_id, getInNeighbors(), getInDegree(), vertexSubset, f, g,
        parallel);
  }

  template <class VS, class F, class G>
  inline void decodeOutNghBreakEarly(uintE vtx_id, VS& vertexSubset, F& f, G& g,
                                     bool parallel = 0) {
    vertex_ops::decodeNghsBreakEarly<symmetric_vertex, W, F, G, VS>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), vertexSubset, f, g,
        parallel);
  }

  template <class F, class G>
  inline void decodeOutNgh(uintE vtx_id, F& f, G& g) {
    vertex_ops::decodeNghs<symmetric_vertex, W, F, G>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), f, g);
  }

  template <class F, class G>
  inline void decodeInNgh(uintE vtx_id, F& f, G& g) {
    vertex_ops::decodeNghs<symmetric_vertex, W, F, G>(
        this, vtx_id, getInNeighbors(), getInDegree(), f, g);
  }

  template <class F, class G, class H>
  inline void decodeOutNghSparse(uintE vtx_id, uintT o, F& f, G& g, H& h, bool parallel=true) {
    vertex_ops::decodeNghsSparse<symmetric_vertex, W, F>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), o, f, g, h);
  }

  template <class F, class G, class H>
  inline void decodeInNghSparse(uintE vtx_id, uintT o, F& f, G& g, H& h, bool parallel=true) {
    vertex_ops::decodeNghsSparse<symmetric_vertex, W, F>(
        this, vtx_id, getInNeighbors(), getInDegree(), o, f, g, h, parallel);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseSeq(uintE vtx_id, uintT o, F& f, G& g) {
    return vertex_ops::decodeNghsSparseSeq<symmetric_vertex, W, F>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeInNghSparseSeq(uintE vtx_id, uintT o, F& f, G& g) {
    return vertex_ops::decodeNghsSparseSeq<symmetric_vertex, W, F>(
        this, vtx_id, getInNeighbors(), getInDegree(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutBlock(uintE vtx_id, uintT o, uintE block_num, F& f,
                               G& g) {
    return vertex_ops::decode_block<W, F>(vtx_id, getOutNeighbors(),
                                          getOutDegree(), o, block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeInBlock(uintE vtx_id, uintT o, uintE block_num, F& f,
                              G& g) {
    return decodeOutBlock(vtx_id, o, block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseBlock(uintE vtx_id, uintT o, uintE block_size,
                                        uintE block_num, F& f, G& g) {
    return vertex_ops::decodeNghsSparseBlock<symmetric_vertex, W, F>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), o, block_size,
        block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeInNghSparseBlock(uintE vtx_id, uintT o, uintE block_size,
                                       uintE block_num, F& f, G& g) {
    return vertex_ops::decodeNghsSparseBlock<symmetric_vertex, W, F>(
        this, vtx_id, getInNeighbors(), getInDegree(), o, block_size, block_num,
        f, g);
  }

  template <class F, class G>
  inline void copyOutNgh(uintE vtx_id, uintT o, F& f, G& g) {
    vertex_ops::copyNghs<symmetric_vertex, W>(this, vtx_id, getOutNeighbors(),
                                             getOutDegree(), o, f, g);
  }

  template <class F, class G>
  inline void copyInNgh(uintE vtx_id, uintT o, F& f, G& g) {
    vertex_ops::copyNghs<symmetric_vertex, W>(this, vtx_id, getInNeighbors(),
                                             getInDegree(), o, f, g);
  }

  inline edge_type get_ith_out_neighbor(uintE vtx_id, size_t i) {
    return vertex_ops::get_ith_neighbor<symmetric_vertex, W>(getOutNeighbors(),
                                                            i);
  }

  inline edge_type get_ith_in_neighbor(uintE vtx_id, size_t i) {
    return vertex_ops::get_ith_neighbor<symmetric_vertex, W>(getInNeighbors(),
                                                            i);
  }

  template <class F>
  inline size_t countOutNgh(uintE vtx_id, F& f, bool parallel = true) {
    return vertex_ops::countNghs<symmetric_vertex, W, F>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), f, parallel);
  }

  template <class F>
  inline size_t countInNgh(uintE vtx_id, F& f, bool parallel = true) {
    return vertex_ops::countNghs<symmetric_vertex, W, F>(
        this, vtx_id, getInNeighbors(), getInDegree(), f, parallel);
  }

  template <class E, class M, class Monoid>
  inline E reduceOutNgh(uintE vtx_id, M& m, Monoid& reduce) {
    return vertex_ops::reduceNghs<symmetric_vertex, W, E, M, Monoid>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), m, reduce);
  }

  template <class E, class M, class Monoid>
  inline E reduceInNgh(uintE vtx_id, M& m, Monoid& reduce) {
    return vertex_ops::reduceNghs<symmetric_vertex, W, E, M, Monoid>(
        this, vtx_id, getInNeighbors(), getInDegree(), m, reduce);
  }

  template <class F>
  inline void mapOutNgh(uintE vtx_id, F& f, bool parallel = true) {
    vertex_ops::mapNghs<symmetric_vertex, W, F>(this, vtx_id, getOutNeighbors(),
                                               getOutDegree(), f, parallel);
  }

  template <class F>
  inline void mapInNgh(uintE vtx_id, F& f, bool parallel = true) {
    vertex_ops::mapNghs<symmetric_vertex, W, F>(this, vtx_id, getInNeighbors(),
                                               getInDegree(), f, parallel);
  }

  template <class P, class O>
  inline void filterOutNgh(uintE vtx_id, P& p, O& out,
                           edge_type* tmp) {
    vertex_ops::filterNghs<symmetric_vertex, W, P, O>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), p, out, tmp);
  }

  template <class P, class O>
  inline void filterInNgh(uintE vtx_id, P& p, O& out,
                          edge_type* tmp) {
    vertex_ops::filterNghs<symmetric_vertex, W, P, O>(
        this, vtx_id, getInNeighbors(), getInDegree(), p, out, tmp);
  }

  template <class P>
  inline size_t packOutNgh(uintE vtx_id, P& p, edge_type* tmp) {
    uintE k = vertex_ops::packNghs<symmetric_vertex, W, P>(
        this, vtx_id, p, getOutNeighbors(), getOutDegree(), tmp);
    return k;
  }

  template <class P>
  inline size_t packInNgh(uintE vtx_id, P& p, edge_type* tmp) {
    uintE k = vertex_ops::packNghs<symmetric_vertex, W, P>(
        this, vtx_id, p, getInNeighbors(), getInDegree(), tmp);
    return k;
  }

  inline size_t calculateOutTemporarySpace() {
    return vertex_ops::calculateTemporarySpace(getOutDegree());
  }

  inline size_t calculateInTemporarySpace() {
    return vertex_ops::calculateTemporarySpace(getInDegree());
  }
  */
};


//Graph is unweighted
struct dyn_arr_sym_graph {
  using edge_type = uintE; //No need for tuple with weight

  pbbslib::dyn_arr<dyn_arr_sym_vertex*> vertices;

  size_t n; //Number of vertices
  size_t m; //Number of edges

  //Called to delete graph
  std::function<void()> deletion_fn;


  //n - initial capacity for vertices
  //Will be resized if necessary
  dyn_arr_sym_graph(std::function<void()> _deletion_fn, size_t n):
    deletion_fn(_deletion_fn),
    n(0), //Graph starts with no vertices or edges
    m(0),
    vertices(pbbslib::dyn_arr<dyn_arr_sym_vertex*>(n)) //What is the best size
  {}

  //Could be nullptr if vertex was deleted
  dyn_arr_sym_vertex* get_vertex(uintE i) {
    return vertices.A[i];
  }

  void del() {
    deletion_fn();
  }
  
  //Vertex starts out without neighbors
  //n - initial capacity of neighbors
  void add_vertex(uintE i, size_t n) {
    vertices.A[i] = new dyn_arr_sym_vertex(i, n);
  }

  //Does not actually remove space because index needed to keep id
  //Assumes that the vertex has no neighbors
  void remove_vertex(uintE i) {
    delete vertices.A[i];
    vertices.A[i] = nullptr;
  }

  //Assumes that the vertices have not been removed
  void add_edge(uintE u, uintE v) {
    get_vertex(u)->add_neighbor(v);
    get_vertex(v)->add_neighbor(u);
  }

  //Assumes that the vertices have not been removed
  void remove_edge(uintE u, uintE v) {
    get_vertex(u)->remove_neighbor(v);
    get_vertex(v)->remove_neighbor(u);
    
  }

  /*
  template <class P>
  uintE pack_neighbors(uintE id, P& p, std::tuple<uintE, W>* tmp) {
    uintE new_degree = get_vertex(id).packOutNgh(id, p, tmp);
    v_data[id].degree = new_degree; // updates the degree //
    return new_degree;
  }
  */

  /*
  // degree must be <= old_degree //
  void decrease_vertex_degree(uintE id, uintE degree) {
    assert(degree <= v_data[id].degree);
    v_data[id].degree = degree;
  }
  */

  /*
  void zero_vertex_degree(uintE id) {
    decrease_vertex_degree(id, 0);
  }
  */

  /*
  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    parallel_for(0, n, [&](size_t i) {
      get_vertex(i).mapOutNgh(i, f, parallel_inner_map);
    }, 1);
  }
  */

  /*
  // F : edge -> edge
  template <class F>
  void alter_edges(F f, bool parallel_inner_map = true) {
    abort(); // unimplemented for CSR //
  }
  */

  /*
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
  */
};
