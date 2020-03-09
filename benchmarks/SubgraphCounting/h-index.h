#include <iostream>
#include <map>
#include <set>
#include <limits>
#include <vector>
#include "ligra/ligra.h"
#include "dynamic_symmetric_graph.h"
#include "ligra/pbbslib/sparse_table.h"
#include "pbbslib/hash_table.h"

// Original Code for HSet
/*
struct HSet {
	int f[1000];
	std::set<int> H;
	std::set<int> P;
	std::set<int> B;
	// A set to maintain H
	std::set<int> C[1000];

	int inc(int e) {
		return change(e,f[e] + 1);
	}

	int dec(int e){
		if(f[e] == 0) return H.size();
		return change(e,f[e] - 1);
	}

	int insert(int e,int x) {
		f[e] = x;
		C[x].insert(e);

		if(x > H.size()) {
			H.insert(e);
			if(B.size() == 0) {
				if(C[H.size() - 1].size() == 0) {
					B = C[H.size()];
					C[H.size()] = std::set<int>();
				}
			}else{
				int cur = *(B.begin());
				B.erase(B.begin());
				H.erase(cur);
				C[H.size()].insert(cur);
			}
		}

		return H.size();
	}
	
	int erase(int e) {
		if(B.find(e) != B.end()) {
			B.erase(e);
		}else{
			C[f[e]].erase(e);
		}
		if(H.find(e) != H.end()) {
			H.erase(e);
			if(C[H.size() + 1].size() == 0) {
				C[H.size() + 1] = B;
				B = std::set<int>();
			}else{
				int cur = *(C[H.size() + 1].begin());
				C[H.size() + 1].erase(cur);
				B.insert(cur);
				H.insert(cur);
			}
		}
		return H.size();
	}
	
	int change(int e,int x) {
		erase(e);
		return insert(e,x);
	}
};
*/

// New Code for HSet
// For hashing sparse table
struct hash_uintE {
  inline size_t operator () (const uintE & a) {return pbbs::hash64_2(a);}
};


struct HSet {
  symmetric_graph<symmetric_vertex, pbbs::empty>* G;
  sparse_table<uintE, pbbs::empty, hash_uintE> H;
  sparse_table<uintE, pbbs::empty, hash_uintE> B;

  sparse_table<uintE, sparse_table<uintE, pbbs::empty, hash_uintE>, hash_uintE> C;

  sparse_table<uintE, pbbs::empty, hash_uintE> empty;

  //n - number of vertices
  //a - threshold for low and high degree vertices
  HSet(size_t n, size_t m) {

    empty = make_sparse_table<uintE, pbbs::empty, hash_uintE>(0, std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());

    auto v_data = pbbs::new_array_no_init<vertex_data>(n);
    auto edges = pbbs::new_array_no_init<std::tuple<uintE, pbbs::empty>>(m);
    G = new symmetric_graph<symmetric_vertex, pbbs::empty>(v_data, n, m, get_deletion_fn(v_data, edges), (std::tuple<uintE, pbbs::empty>*) edges);


    H = make_sparse_table<uintE, pbbs::empty, hash_uintE>
      (n, std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());

    B = make_sparse_table<uintE, pbbs::empty, hash_uintE>
      (n, std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());

    C = make_sparse_table<uintE, sparse_table<uintE, pbbs::empty, hash_uintE>, hash_uintE>
      (n, std::make_tuple(UINT_E_MAX, empty), hash_uintE());
  }

  uintE insert(uintE v) {

    //ADD TO GRAPH!!
    auto deg = G->get_vertex(v).getOutDegree();
    addToC(v);

    if(deg > H.sizeOf()) {

      H.insert(std::make_tuple(v, pbbs::empty()));

      if(B.sizeOf() != 0) {
        auto y = *(B.begin());
        B.erase(y);
        H.erase(y);
        addToC(y);
      }

      else {
        // If C has entry for |H|, set that to B. Otherwise B is empty
        if (C.contains(H.sizeOf())) {

          B = C.find(H.sizeOf(), empty);
          C.erase(H.sizeOf());
        }
        else {
          B.clear();
        }
      }

    }

    return H.sizeOf();
  }

   
  uintE erase(uintE v) {

    // Remove from Graph!!
    auto deg = G->get_vertex(v).getOutDegree();
    if (B.contains(v)) {
      B.erase(v);
    }
    else {
      C.find(deg, empty).erase(v);
    }

    
    if (H.contains(v)) {
      auto h = H.sizeOf(); //h is |H| before removing v
      H.erase(v);

      if (C.find(h, empty).sizeOf() != 0) { //Has entry for C[h]
        auto y = *C.find(h, empty).begin();
        removeFromC(y);
        B.insert(std::make_tuple(y, pbbs::empty()));
        H.insert(std::make_tuple(y, pbbs::empty()));
      }
      else {
        C.insert(std::make_tuple(h, B));
        B.clear();
      }

    }

    return H.sizeOf();
  }
	
  uintE change(uintE v, uintE x) {
    erase(v);
    return insert(v);
  }

  void addToC(uintE v) {

    uintE deg = G->get_vertex(v).getOutDegree();

    auto vertices = C.find(deg, empty); //Vertices of degree "deg"
	    
    if (vertices.sizeOf() == 0) { //If entry for deg is empty, create new sparse_table entry
        C.insert(std::make_tuple(deg, make_sparse_table<uintE, pbbs::empty, hash_uintE>(G->n, std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE())));
    }
    C.find(deg, empty).insert(std::make_tuple(deg, pbbs::empty()));

  }

  void removeFromC(uintE v) { // Assumes v is in C
    
    uintE deg = G->get_vertex(v).getOutDegree();

    auto vertices = C.find(deg, empty);
    if (vertices.sizeOf() == 1) { // Deletes entire entry if it is empty (or will be)
      C.erase(deg);
    }
    else {
      vertices.erase(v);
    }
  }

};


struct graph {

  public:
    HSet* h;
    // TODO: don't use constants like these; should be malloc-ed
    // using sequence would probably be easier than using malloc directly
    // Also, edges should come from the passed-in Graph template object
    std::vector<int> edges[1000];
    std::vector<int> hedges[100][100];

	public:
		HSet h;
		std::vector<int> edges[1000];
		std::vector<int> hedges[100][100];

    int triangles = 0;

    graph() {
      h = new HSet(1000, 100);
    }

    void insertV(uintE v) {
      h->insert(v);
    }

    int connect(uintE v, uintE u) {
      // Add u and v to H
      //h.inc(v);
      //h.inc(u);

      triangles += hedges[u][v].size();
      for(int i = 0;i < edges[u].size();i++) {
        // Add new wedge v -- u -- i
        hedges[v][edges[u][i]].push_back(u);
        hedges[edges[u][i]][v].push_back(u);
      }
      for(int i = 0;i < edges[v].size();i++) {
        // Add new wedge u -- v -- i
        hedges[u][edges[v][i]].push_back(v);
        hedges[edges[v][i]][u].push_back(v);
      }
      edges[v].push_back(u);
      edges[u].push_back(v);

      return triangles;
    }
};
		int connect(int v,int u) {
			h.inc(v);
			h.inc(u);

			triangles += hedges[u][v].size();
			for(int i = 0;i < edges[u].size();i++) {
				hedges[v][edges[u][i]].push_back(u);
				hedges[edges[u][i]][v].push_back(u);
			}
			for(int i = 0;i < edges[v].size();i++) {
				hedges[u][edges[v][i]].push_back(v);
				hedges[edges[v][i]][u].push_back(v);
			}
			edges[v].push_back(u);
			edges[u].push_back(v);

			return this -> triangles;
		}
};
