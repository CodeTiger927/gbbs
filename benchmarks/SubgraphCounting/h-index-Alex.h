
#include <iostream>
#include <limits>
#include <vector>
#include <math.h>
#include "ligra/ligra.h"
#include "dynamic_symmetric_graph.h"
#include "ligra/pbbslib/sparse_table.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "pbbslib/sequence.h"

const size_t SIZE_OF_GRAPH = 10;

struct HSetAlex {
	HSetAlex() {
		init();
	}

	uintE n;
	uintE HIndex;
	sparse_table<uintE,bool,hash_uintE> B;
	uintE BSize;
	pbbslib::dyn_arr<sparse_table<uintE,bool,hash_uintE>> C;
	pbbslib::dyn_arr<uintE> cN;
	pbbslib::dyn_arr<uintE> vertices;

	// This can potentially be O(N), but I think normally it would be smaller than O(h) in actual practice
	sequence<uintE> allH() {
		auto all = make_sparse_table<uintE,bool,hash_uintE>(2 * HIndex,std::make_tuple(UINT_E_MAX,false),hash_uintE());
		par_for(HIndex + 1,n,[&](size_t i) {
			if(cN.A[i] != 0) {
				auto entries = C.A[i].entries();
				par_for(0,entries.size(),[&](size_t j) {
					all.insert(std::make_tuple(std::get<0>(entries[j]),true));
				});
			}
		});

		auto b = B.entries();
		par_for(0,b.size(),[&](size_t i) {
			all.insert(std::make_tuple(std::get<0>(b[i]),true));
		});

		auto f = all.entries();

		return sequence<uintE>(f.size(),[&](size_t i) {return std::get<0>(f[i]);});
	}

	bool inH(uintE v) {
		if(vertices.A[v] > HIndex) return true;
		if(vertices.A[v] == HIndex) return B.find(v,false);
		return false;
	}

	void init() {
		n = 0;
		C = pbbslib::dyn_arr<sparse_table<uintE,bool,hash_uintE>>(0);
		cN = pbbslib::dyn_arr<uintE>(0);
		vertices = pbbslib::dyn_arr<uintE>(0);
		HIndex = 0;
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
      uintE* nA = pbbslib::new_array_no_init<uintE>(nC);
      par_for(0,n,1,[&](size_t i){nA[i] = vertices.A[i];});
      pbbslib::free_array(vertices.A);
      vertices.A = nA;
      vertices.capacity = nC;

      sparse_table<uintE,bool,hash_uintE>* nAA = pbbslib::new_array_no_init<sparse_table<uintE,bool,hash_uintE>>(nC);
      par_for(0,n,1,[&](size_t i){nAA[i] = C.A[i];});
      pbbslib::free_array(C.A);
      C.A = nAA;
      C.capacity = nC;

      uintE* nAAA = pbbslib::new_array_no_init<uintE>(nC);
      par_for(0,n,1,[&](size_t i){nAAA[i] = cN.A[i];});
      pbbslib::free_array(cN.A);
      cN.A = nAAA;
      cN.capacity = nC;

      n = amount;
   }

	void resizeC(uintE v,uintE amount) {
		amount = std::max(amount,(uintE)SIZE_OF_GRAPH);
		if(amount != 0 && ceil(log2(amount)) == ceil(log2(C.A[v].m))) return;
		auto entries = C.A[v].entries();

		C.A[v] = make_sparse_table<uintE,bool,hash_uintE>(amount,std::make_tuple(UINT_E_MAX,false),hash_uintE());
		par_for(0,entries.size(),1,[&](size_t i) {C.A[v].insert(entries[i]);});
	}

	void remove(sequence<uintE> s) {
		n -= s.size();
		sequence<uintE> all = merge_sort(s,[&](uintE a,uintE b) {return vertices.A[a] < vertices.A[b];});
		long long curN = HIndex - BSize + cN.A[HIndex];
    	par_for(0,s.size(),[&](size_t i) {
			C.A[vertices.A[s[i]]].erase(s[i]);
		});

    	auto start = make_sparse_table<uintE,uintE,hash_uintE>(s.size(),std::make_tuple(UINT_E_MAX,UINT_E_MAX),hash_uintE());
		par_for(0,all.size(),1,[&](size_t i) {
			if(i != 0) {
				if(vertices.A[all[i]] != vertices.A[all[i - 1]]) {
					start.insert(std::make_tuple(vertices.A[all[i]],i));
					cN.A[vertices.A[all[i - 1]]] -= i;
				}
			}
		});


		start.insert(std::make_tuple(vertices.A[all[0]],0));
		cN.A[vertices.A[all[all.size() - 1]]] -= all.size();
		auto entries = start.entries();


		par_for(0,entries.size(),[&](size_t i) {
			uintE cur = std::get<0>(entries[i]);
			cN.A[cur] += std::get<1>(entries[i]);
			resizeC(cur,cN.A[cur]);
		});


		curN -= filter(s,[&](uintE i){return vertices.A[i] >= HIndex;}).size();
      	par_for(0,s.size(),[&](size_t i) {
         	vertices.A[s[i]] = 0;
      	});

		sequence<uintE> cNs = sequence<uintE>(s.size() + 2,[&](size_t i) {
			long long cur = HIndex - i;
			if(cur < 0) return (uintE)0;
			return cN.A[cur];
		});



		pbbslib::scan_add_inplace(cNs);
		sequence<uintE> worksOrNo = sequence<uintE>(cNs.size() + 1,[&](size_t i) {
			if(curN + cNs[i] >= HIndex - i) {
				return (uintE)(HIndex - i);
			}else{
				return (uintE)0;
			}
		});

		uintE NHIndex = pbbs::reduce(worksOrNo,pbbs::maxm<uintE>());
		if(NHIndex == HIndex) return;
		auto allH = C.A[HIndex].entries();
		BSize = NHIndex - (curN - cNs[HIndex - NHIndex + 1]);
		B = make_sparse_table(std::max(BSize,(uintE)SIZE_OF_GRAPH),std::make_tuple(UINT_E_MAX,false),hash_uintE());
		par_for(0,BSize,[&](size_t i) {
			B.insert(allH[i]);
		});

		HIndex = NHIndex;
	}

	void insert(sequence<std::pair<uintE,uintE>> s) {
		n += s.size();
    	sequence<std::pair<uintE,uintE>> all = merge_sort(s,[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a.second > b.second;});
		long long curN = HIndex - BSize + cN.A[HIndex];
    	auto start = make_sparse_table<uintE,uintE,hash_uintE>(s.size(),std::make_tuple(UINT_E_MAX,UINT_E_MAX),hash_uintE());
		par_for(0,s.size(),1,[&](size_t i) {
			vertices.A[s[i].first] = s[i].second;
			if(i != 0) {
				if(all[i].second != all[i - 1].second) {
					start.insert(std::make_tuple(all[i].second,i));
					cN.A[all[i - 1].second] += i;
				}
			}
		});

		start.insert(std::make_tuple(all[0].second,0));
		cN.A[all[all.size() - 1].second] += all.size();
		auto entries = start.entries();

		par_for(0,entries.size(),[&](size_t i) {
			uintE cur = std::get<0>(entries[i]);
			cN.A[std::get<0>(entries[i])] -= std::get<1>(entries[i]);
			resizeC(std::get<0>(entries[i]),cN.A[cur]);
		});

		par_for(0,s.size(),[&](size_t i) {
			C.A[s[i].second].insert(std::make_tuple(s[i].first,true));
		});


		curN += filter(s,[&](std::pair<uintE,uintE> i){return i.second >= HIndex;}).size();
		sequence<uintE> cNs = sequence<uintE>(s.size() + 2,[&](size_t i) {
			uintE cur = i + HIndex;
			if(cur >= SIZE_OF_GRAPH) return (uintE)0;
			return cN.A[cur];
		});


		pbbslib::scan_add_inplace(cNs);

		sequence<uintE> worksOrNo = sequence<uintE>(cNs.size(),[&](size_t i) {
			if(curN - cNs[i] >= HIndex + i) {
				return (uintE)(HIndex + i);
			}else{
				return HIndex;
			}
		});


		uintE NHIndex = pbbs::reduce(worksOrNo,pbbs::maxm<uintE>());
		if(NHIndex == HIndex) return;
		auto allH = C.A[NHIndex].entries();
		BSize = NHIndex - (curN - cNs[NHIndex - HIndex + 1]);
		B = make_sparse_table(2 * BSize,std::make_tuple(UINT_E_MAX,false),hash_uintE());

		par_for(0,BSize,[&](size_t i) {
			B.insert(allH[i]);
		});
		HIndex = NHIndex;
	}

   void modify(sequence<std::pair<uintE,uintE>> s) {
      sequence<uintE> tR = sequence<uintE>(s.size(),[&](size_t i) {return s[i].first;});
      remove(tR);
      insert(s);
   }
};