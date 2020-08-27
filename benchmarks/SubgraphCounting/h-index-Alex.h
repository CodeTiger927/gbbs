
#include <iostream>
#include <limits>
#include <vector>
#include <math.h>
#include "ligra/ligra.h"
#include "dynamic_symmetric_graph.h"
#include "ligra/pbbslib/sparse_table.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "pbbslib/sequence.h"

size_t SIZE_OF_GRAPH = 1005;

class HSet {

  public:
    dynamic_symmetric_graph<dynamic_symmetric_vertex, uintE>* G;
    size_t hindex;

    HSet(dynamic_symmetric_graph<dynamic_symmetric_vertex, uintE>* _G) {
      G = _G;
      hindex = 0;
    }

    virtual pbbs::sequence<uintE> getH() = 0;
    virtual bool contains(uintE target) = 0;

    virtual uintE insertVertices(pbbs::sequence<uintE> vertices) = 0;
    virtual uintE eraseVertices(pbbs::sequence<uintE> vertices) = 0;
    virtual uintE insertEdges(pbbs::sequence<std::pair<uintE, uintE>> edges) = 0;
    virtual uintE eraseEdges(pbbs::sequence<std::pair<uintE, uintE>> edges) = 0;
};


class HSetAlex : public HSet {
public:
	HSetAlex(dynamic_symmetric_graph<dynamic_symmetric_vertex, uintE>* graph): HSet(graph) {
		init();
	}
	uintE n;
	sparse_table<uintE,bool,hash_uintE> B;
	uintE BSize;
	pbbslib::dyn_arr<sparse_table<uintE,bool,hash_uintE>> C;
	pbbslib::dyn_arr<uintE> cN;
	pbbslib::dyn_arr<uintE> cStored;

	// This can potentially be O(N), but I think normally it would be smaller than O(h) in actual practice
	sequence<uintE> allH() {
		// Make this O(H) by making a sparse table on all the degrees with at least 1. 
		auto all = make_sparse_table<uintE,bool,hash_uintE>(2 * hindex + 1,std::make_tuple(UINT_E_MAX,false),hash_uintE());
		
		par_for(hindex + 1,550,[&](size_t i) {
			if(cN.A[i] != 0) {
				auto entries = C.A[i].entries();
				par_for(0,entries.size(),100,[&](size_t j) {
					if(std::get<1>(entries[j])) all.insert(std::make_tuple(std::get<0>(entries[j]),true));
				});
			}
		});

		auto b = B.entries();

		par_for(0,b.size(),[&](size_t i) {
			if(std::get<1>(b[i])) all.insert(std::make_tuple(std::get<0>(b[i]),true));
		});

		auto f = all.entries();

		return sequence<uintE>(f.size(),[&](size_t i) {return std::get<0>(f[i]);});
	}
	sequence<uintE> getH() {
		return allH();
	}

	bool inH(uintE v) {
		if(G -> v_data.A[v].degree > hindex) return true;
		if(G -> v_data.A[v].degree == hindex) return B.find(v,false);
		return false;
	}

	bool contains(uintE v) {
		return inH(v);
	}

	void init() {
		n = 0;
		C = pbbslib::dyn_arr<sparse_table<uintE,bool,hash_uintE>>(0);
		cN = pbbslib::dyn_arr<uintE>(0);
		cStored = pbbslib::dyn_arr<uintE>(0);

		hindex = 0;
		BSize = 0;
		B = make_sparse_table<uintE,bool,hash_uintE>(SIZE_OF_GRAPH,std::make_tuple(UINT_E_MAX,false),hash_uintE());
	}

   void resizeV(size_t amount) {
	  amount = std::max(amount,INIT_DYN_GRAPH_EDGE_SIZE);
	  size_t cur = n;
	  if(amount != 0 && ceil(log2(amount)) == ceil(log2(cur))) {
		 return;
	  }

	  size_t nC = 1 << ((size_t)(ceil(log2(amount))));

	  sparse_table<uintE,bool,hash_uintE>* nAA = pbbslib::new_array_no_init<sparse_table<uintE,bool,hash_uintE>>(nC);
	  par_for(0,n,1,[&](size_t i){nAA[i] = C.A[i];});
	  par_for(n,nC,1,[&](size_t i){nAA[i] = make_sparse_table<uintE,bool,hash_uintE>(1,std::make_tuple(UINT_E_MAX,false),hash_uintE());});
	  
	  pbbslib::free_array(C.A);
	  C.A = nAA;
	  C.capacity = nC;


	  uintE* nAAA = pbbslib::new_array_no_init<uintE>(nC);
	  par_for(0,n,1,[&](size_t i){nAAA[i] = cN.A[i];});
	  par_for(n,nC,1,[&](size_t i){nAAA[i] = 0;});
	  pbbslib::free_array(cN.A);
	  cN.A = nAAA;
	  cN.capacity = nC;

	  uintE* nAAAA = pbbslib::new_array_no_init<uintE>(nC);
	  par_for(0,n,1,[&](size_t i){nAAAA[i] = cStored.A[i];});
	  par_for(n,nC,1,[&](size_t i){nAAAA[i] = 0;});
	  pbbslib::free_array(cStored.A);
	  cStored.A = nAAAA;
	  cStored.capacity = nC;

	  n = amount;
	  SIZE_OF_GRAPH = amount;
   }

	void resizeC(uintE v) {
		uintE amount = std::max(cStored.A[v],(uintE)SIZE_OF_GRAPH);
		if((amount << 1) <= C.A[v].m && (amount << 2) >= C.A[v].m) return; 
		auto entries = C.A[v].entries();
		C.A[v] = make_sparse_table<uintE,bool,hash_uintE>(2 * amount,std::make_tuple(UINT_E_MAX,false),hash_uintE());
		par_for(0,entries.size(),1,[&](size_t i) {if(std::get<1>(entries[i])) C.A[v].insert(entries[i]);});
	}

	void remove(sequence<uintE> s) {
		n -= s.size();

		sequence<uintE> all = merge_sort(s,[&](uintE a,uintE b) {return G -> v_data.A[a].degree < G -> v_data.A[b].degree;});
		
		long long curN = hindex - BSize + cN.A[hindex];
		par_for(0,s.size(),[&](size_t i) {
			C.A[G -> v_data.A[s[i]].degree].change(s[i],false);
		});

		auto start = make_sparse_table<uintE,uintE,hash_uintE>(2 * s.size() + 1,std::make_tuple(UINT_E_MAX,UINT_E_MAX),hash_uintE());
		par_for(0,all.size(),1,[&](size_t i) {
			if(i != 0) {
				if(G -> v_data.A[all[i]].degree != G -> v_data.A[all[i - 1]].degree) {
					start.insert(std::make_tuple(G -> v_data.A[all[i]].degree,i));
					cN.A[G -> v_data.A[all[i - 1]].degree] -= i;
				}
			}
		});
		start.insert(std::make_tuple(G -> v_data.A[all[0]].degree,0));
		cN.A[G -> v_data.A[all[all.size() - 1]].degree] -= all.size();
		auto entries = start.entries();


		par_for(0,entries.size(),[&](size_t i) {
			uintE cur = std::get<0>(entries[i]);
			cN.A[cur] += std::get<1>(entries[i]);
			resizeC(cur);
		});

		curN -= filter(s,[&](uintE i){return  G -> v_data.A[i].degree >= hindex;}).size();

		sequence<uintE> cNs = sequence<uintE>(s.size() + 2,[&](size_t i) {
			long long cur = hindex - i - 1;
			if(cur < 0) return (uintE)0;
			return cN.A[cur];
		});

		pbbslib::scan_add_inplace(cNs);
		sequence<long long> worksOrNo = sequence<long long>(cNs.size() + 1,[&](size_t i) {
			if(((long long)curN + (long long)cNs[i]) >= ((long long)hindex - (long long)i)) {
				return ((long long)hindex - (long long)i);
			}else{
				return (long long)0;
			}
		});

		uintE Nhindex = pbbs::reduce(worksOrNo,pbbs::maxm<long long>());
		BSize = cN.A[Nhindex] - (curN + cNs[hindex - Nhindex] - Nhindex);
		auto allH = C.A[Nhindex].entries();		

		B = make_sparse_table(2 * std::max(BSize,(uintE)SIZE_OF_GRAPH) + 1,std::make_tuple(UINT_E_MAX,false),hash_uintE());

		par_for(0,BSize,[&](size_t i) {
			B.insert(allH[i]);
		});


		hindex = Nhindex;
	}

	void insert(sequence<std::pair<uintE,uintE>> s) {

		n += s.size();
		sequence<std::pair<uintE,uintE>> all = merge_sort(s,[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a.second > b.second;});
		long long curN = hindex - BSize + cN.A[hindex];

		auto start = make_sparse_table<uintE,uintE,hash_uintE>(s.size(),std::make_tuple(UINT_E_MAX,UINT_E_MAX),hash_uintE());
		par_for(0,s.size(),1,[&](size_t i) {
			if(i != 0) {
				if(all[i].second != all[i - 1].second) {
					start.insert(std::make_tuple(all[i].second,i));
					cN.A[all[i - 1].second] += i;
					cStored.A[all[i - 1].second] += i;
				}
			}
		});

		start.insert(std::make_tuple(all[0].second,0));
		cN.A[all[all.size() - 1].second] += all.size();
		cStored.A[all[all.size() - 1].second] += all.size();
		auto entries = start.entries();

		par_for(0,entries.size(),[&](size_t i) {
			uintE cur = std::get<0>(entries[i]);
			cN.A[std::get<0>(entries[i])] -= std::get<1>(entries[i]);
			cStored.A[std::get<0>(entries[i])] -= std::get<1>(entries[i]);
			resizeC(std::get<0>(entries[i]));
		});

		par_for(0,s.size(),[&](size_t i) {
			C.A[s[i].second].insert(std::make_tuple(s[i].first,true));
		});

		curN += filter(s,[&](std::pair<uintE,uintE> i){return i.second >= hindex;}).size();
		sequence<uintE> cNs = sequence<uintE>(s.size() + 2,[&](size_t i) {
			uintE cur = i + hindex;
			if(cur >= SIZE_OF_GRAPH) return (uintE)0;
			return cN.A[cur];
		});

		pbbslib::scan_add_inplace(cNs);


		sequence<uintE> worksOrNo = sequence<uintE>(cNs.size(),[&](size_t i) {
			if((long long)curN - (long long)cNs[i] >= (long long)hindex + (long long)i) {
				return (uintE)(hindex + i);
			}else{
				return (uintE)hindex;
			}
		});


		uintE Nhindex = pbbs::reduce(worksOrNo,pbbs::maxm<uintE>());

		auto allH = C.A[Nhindex].entries();
		BSize = Nhindex - (curN - cNs[Nhindex - hindex + 1]);
		B = make_sparse_table(2 * BSize + 1,std::make_tuple(UINT_E_MAX,false),hash_uintE());

		par_for(0,BSize,[&](size_t i) {
			B.insert(allH[i]);
		});
		hindex = Nhindex;
	}

   void modify(sequence<std::pair<uintE,uintE>> s) {
		sequence<uintE> tR = sequence<uintE>(s.size(),[&](size_t i) {return s[i].first;});
		remove(tR);
		insert(s);
   }

	uintE insertVertices(pbbs::sequence<uintE> vertices) {
		G -> batchAddVertices(vertices);
	}

	uintE eraseVertices(pbbs::sequence<uintE> vertices) {
		G -> batchRemoveVertices(vertices);
	}

	uintE insertEdges(pbbs::sequence<std::pair<uintE, uintE>> edges) {
		auto aN = merge_sort(sequence<uintE>(2 * edges.size(),[&](size_t i) {
			if(i % 2) {
				return edges[i / 2].first;
			}else{
				return edges[i / 2].second;
			}}),[&](uintE a,uintE b){return a > b;});

		auto allNodes = sequence<uintE>(2 * edges.size(),[&](size_t i) {
			if(i == 0) {
				return aN[i];
			}else{
				if(aN[i] != aN[i - 1]) {
					return aN[i];
				}else{
					return UINT_E_MAX;
				}
			}
  		});

		allNodes = filter(allNodes,[&](uintE i) {return (i != UINT_E_MAX);});
		// for(int i = 0;i < allNodes.size();++i) {
		// 	cout << allNodes[i] << " " << G -> v_data.A[allNodes[i]].degree << endl;
		// }
		remove(allNodes);
		G -> batchAddEdges(edges);
		insert(sequence<std::pair<uintE,uintE>>(allNodes.size(),[&](size_t i) {return std::make_pair(allNodes[i],G -> v_data.A[allNodes[i]].degree);}));
	}

	uintE eraseEdges(pbbs::sequence<std::pair<uintE, uintE>> edges) {
		auto aN = merge_sort(sequence<uintE>(2 * edges.size(),[&](size_t i) {
			if(i % 2) {
				return edges[i / 2].first;
			}else{
				return edges[i / 2].second;
		}}),[&](uintE a,uintE b){return a > b;});

		auto allNodes = sequence<uintE>(2 * edges.size(),[&](size_t i) {
			if(i == 0) {
				return aN[i];
			}else{
				if(aN[i] != aN[i - 1]) {
					return aN[i];
				}else{
					return UINT_E_MAX;
				}
			}
		});

		allNodes = filter(allNodes,[&](uintE i) {return (i != UINT_E_MAX);});

		remove(allNodes);
		
		G -> batchRemoveEdges(edges);

		insert(sequence<std::pair<uintE,uintE>>(allNodes.size(),[&](size_t i) {return std::make_pair(allNodes[i],G -> v_data.A[allNodes[i]].degree);}));
    }
};