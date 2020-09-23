
#include <iostream>
#include <limits>
#include <vector>
#include <math.h>
#include "ligra/ligra.h"
#include "dynamic_symmetric_graph.h"
#include "ligra/pbbslib/sparse_table.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "pbbslib/sequence.h"
#include "hindex.h"

class HSetSparseTable : public HSet {
public:
	HSetSparseTable(dynamic_symmetric_graph<dynamic_symmetric_vertex, pbbs::empty>* graph, uintE size): HSet(graph) {
		initExtra(size);
	}
	uintE n;
	sparse_table<uintE,bool,hash_uintE> B;
	sparse_table<uintE,bool,hash_uintE> eDegs;
	uintE BSize;
	pbbslib::dyn_arr<sparse_table<uintE,bool,hash_uintE>> C;
	pbbslib::dyn_arr<uintE> cN;
	pbbslib::dyn_arr<uintE> cStored;

	// This can potentially be O(N), but I think normally it would be smaller than O(h) in actual practice
	sequence<uintE> allH() {
		// Make this O(H) by making a sparse table on all the degrees with at least 1. 
		auto all = make_sparse_table<uintE,bool,hash_uintE>(2 * hindex + 1,std::make_tuple(UINT_E_MAX,false),hash_uintE());
		// par_for(hindex + 1,550,[&](size_t i) {
		// 	if(cN.A[i] != 0) {
		// 		auto entries = C.A[i].entries();
		// 		par_for(0,entries.size(),100,[&](size_t j) {
		// 			if(std::get<1>(entries[j])) all.insert(std::make_tuple(std::get<0>(entries[j]),true));
		// 		});
		// 	}
		// });
		auto allDegs = eDegs.entries();
		par_for(0,allDegs.size(),[&](size_t i) {
			if(!std::get<1>(allDegs[i]) || std::get<0>(allDegs[i]) <= hindex) return;
			auto entries = C.A[std::get<0>(allDegs[i])].entries();
			par_for(0,entries.size(),1,[&](size_t j) {
				if(std::get<1>(entries[j])) all.insert(std::make_tuple(std::get<0>(entries[j]),true));
			});
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
		B = make_sparse_table<uintE,bool,hash_uintE>(INIT_C_SIZE,std::make_tuple(UINT_E_MAX,false),hash_uintE());
		eDegs = make_sparse_table<uintE,bool,hash_uintE>(INIT_C_SIZE,std::make_tuple(UINT_E_MAX,false),hash_uintE());
	}

	void initExtra(uintE n) {
		resizeV(n);
		insert(sequence<std::pair<uintE,uintE>>(n,[&](size_t i) {return std::make_pair(i,0);}));
	}

	void resizeDegs(size_t size) {
		auto entries = eDegs.entries();
		eDegs = make_sparse_table<uintE,bool,hash_uintE>(size,std::make_tuple(UINT_E_MAX,false),hash_uintE());
		par_for(0,entries.size(),[&](size_t i) {if(std::get<1>(entries[i])) eDegs.insert(entries[i]);});
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
	  par_for(n,nC,1,[&](size_t i){nAA[i] = make_sparse_table<uintE,bool,hash_uintE>(INIT_C_SIZE,std::make_tuple(UINT_E_MAX,false),hash_uintE());});
	  
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
   }

	void resizeC(uintE v) {
		uintE amount = std::max(cStored.A[v],(uintE)INIT_DEG_SIZE);
		if(amount << 1 <= C.A[v].m && ((std::max(cN.A[v],(uintE)INIT_DEG_SIZE)) << 2) >= C.A[v].m) return; 
		auto entries = C.A[v].entries();
		// Something weird is going on here. Not sure what
		C.A[v] = make_sparse_table<uintE,bool,hash_uintE>(2 * (std::max(cN.A[v],(uintE)INIT_C_SIZE)),std::make_tuple(UINT_E_MAX,false),hash_uintE());
		par_for(0,entries.size(),1,[&](size_t i) {if(std::get<1>(entries[i])) C.A[v].insert(entries[i]);});
		cStored.A[v] = cN.A[v];
	}

	void remove(sequence<uintE> s) {
		n -= s.size();

		sequence<uintE> all = merge_sort(s,[&](uintE a,uintE b) {return G -> v_data.A[a].degree < G -> v_data.A[b].degree;});
		
		long long curN = hindex - BSize + cN.A[hindex];
		par_for(0,s.size(),[&](size_t i) {
			cout << G -> v_data.A[s[i]].degree << " " << C.A[G -> v_data.A[s[i]].degree].m << " " << s[i] << endl;
			C.A[G -> v_data.A[s[i]].degree].change(s[i],false);
		});

		par_for(0,all.size(),1,[&](size_t i) {
			if(i != 0) {
				if(G -> v_data.A[all[i]].degree != G -> v_data.A[all[i - 1]].degree) {
					cN.A[G -> v_data.A[all[i - 1]].degree] -= i;
				}
			}
		});

		cN.A[G -> v_data.A[all[all.size() - 1]].degree] -= all.size();

		par_for(0,all.size(),1,[&](size_t i) {
			if(i != 0) {
				if(G -> v_data.A[all[i]].degree != G -> v_data.A[all[i - 1]].degree) {
					cN.A[G -> v_data.A[all[i]].degree] += i;
					resizeC(G -> v_data.A[all[i]].degree);
				}
			}
		});

		resizeC(G -> v_data.A[all[0]].degree);

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

		B = make_sparse_table(2 * std::max(BSize,(uintE)INIT_C_SIZE) + 1,std::make_tuple(UINT_E_MAX,false),hash_uintE());

		par_for(0,BSize,[&](size_t i) {
			B.insert(allH[i]);
		});

		hindex = Nhindex;

		par_for(0,s.size(),[&](size_t i) {
			if(cN.A[G -> v_data.A[s[i]].degree] == 0) eDegs.change(G -> v_data.A[s[i]].degree,false);
		});
		resizeDegs(2 * hindex);

	}

	void insert(sequence<std::pair<uintE,uintE>> s) {

		n += s.size();
		sequence<std::pair<uintE,uintE>> all = merge_sort(s,[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a.second > b.second;});
		long long curN = hindex - BSize + cN.A[hindex];

		par_for(0,s.size(),1,[&](size_t i) {
			if(i != 0) {
				if(all[i].second != all[i - 1].second) {
					cN.A[all[i - 1].second] += i;
					cStored.A[all[i - 1].second] += i;
				}
			}
		});

		cN.A[all[all.size() - 1].second] += all.size();
		cStored.A[all[all.size() - 1].second] += all.size();

		par_for(0,s.size(),1,[&](size_t i) {
			if(i != 0) {
				if(all[i].second != all[i - 1].second) {
					cN.A[all[i].second] -= i;
					cStored.A[all[i].second] -= i;
					resizeC(all[i].second);
				}
			}
		});

		resizeC(all[0].second);

		par_for(0,s.size(),[&](size_t i) {
			C.A[s[i].second].insert(std::make_tuple(s[i].first,true));
		});

		curN += filter(s,[&](std::pair<uintE,uintE> i){return i.second >= hindex;}).size();
		sequence<uintE> cNs = sequence<uintE>(s.size() + 2,[&](size_t i) {
			uintE cur = i + hindex;
			if(cur >= n) return (uintE)0;
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

		resizeDegs(2 * hindex);
		par_for(0,s.size(),[&](size_t i) {
			eDegs.insert(std::make_tuple(G -> v_data.A[s[i].first].degree,true));
			eDegs.change(G -> v_data.A[s[i].first].degree,true);
			eDegs.insert(std::make_tuple(G -> v_data.A[s[i].second].degree,true));
			eDegs.change(G -> v_data.A[s[i].second].degree,true);
		});
	}

   void modify(sequence<std::pair<uintE,uintE>> s) {
		sequence<uintE> tR = sequence<uintE>(s.size(),[&](size_t i) {return s[i].first;});
		remove(tR);
		insert(s);
   }

	uintE insertVertices(pbbs::sequence<uintE> vertices) {
		G -> batchAddVertices(vertices);
		return this->hindex;
	}

	uintE eraseVertices(pbbs::sequence<uintE> vertices) {
		G -> batchRemoveVertices(vertices);
		return this->hindex;
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

		return this->hindex;
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

		return this->hindex;
    }
};
