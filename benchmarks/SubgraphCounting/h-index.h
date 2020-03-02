#include <iostream>
#include <map>
#include <set>
#include <limits>
#include <vector>
#include "ligra/ligra.h"
#include "ligra/pbbslib/sparse_table.h"
#include "pbbslib/hash_table.h"

struct HSet {
	// TODO: don't use a constant like 1000; it should be passed in
	// Here, f actually is just the degree, so we don't need it at all; just pass
	// in and save a copy of the graph, so HSet should be templatized
	int f[1000];
	// TODO: replace with GBBS
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


// TODO: don't call this graph -- it's confusing with the Graph templates
// Separate out the way we store the graph from the triangle-counting code
// Also, this isn't quite what we need? It seems like there's no separation of
// H (H-index) and vertices not in H.
struct graph {
	public:
		HSet h;
		// TODO: don't use constants like these; should be malloc-ed
		// using sequence would probably be easier than using malloc directly
		// Also, edges should come from the passed-in Graph template object
		std::vector<int> edges[1000];
		std::vector<int> hedges[100][100];

		int triangles = 0;

		graph() {

		}

		void insertV(int v) {
			h.insert(v,0);
		}

		int connect(int v,int u) {
			// Add u and v to H
			h.inc(v);
			h.inc(u);

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