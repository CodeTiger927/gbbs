
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

	void init() {
		n = 0;
		C = pbbslib::dyn_arr<sparse_table<uintE,bool,hash_uintE>>(SIZE_OF_GRAPH);
		cN = pbbslib::dyn_arr<uintE>(SIZE_OF_GRAPH);
		HIndex = 0;
		BSize = 0;
		B = make_sparse_table<uintE,bool,hash_uintE>(SIZE_OF_GRAPH,std::make_tuple(UINT_E_MAX,false),hash_uintE());
	}

	void resizeC(uintE v,uintE amount) {
		if(amount != 0 && ceil(log2(amount)) == ceil(log2(C.A[v].m))) return;
		auto entries = C.A[v].entries();
		C.A[v] = make_sparse_table<uintE,bool,hash_uintE>(amount,std::make_tuple(UINT_E_MAX,false),hash_uintE());
		par_for(0,entries.size(),1,[&](size_t i) {C.A[v].insert(entries[i]);});
	}

	void insert(sequence<std::pair<uintE,uintE>> s) {
    	sequence<std::pair<uintE,uintE>> all = merge_sort(s,[&](std::pair<uintE,uintE> a,std::pair<uintE,uintE> b) {return a.second > b.second;});
		long long curN = HIndex - BSize + cN.A[HIndex];
    	auto start = make_sparse_table<uintE,uintE,hash_uintE>(s.size(),std::make_tuple(UINT_E_MAX,UINT_E_MAX),hash_uintE());
		par_for(0,s.size(),1,[&](size_t i) {
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
			cN.A[cur] -= std::get<1>(entries[i]);
			resizeC(cur,cN.A[cur]);
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
			if(curN - cNs[i + 1] >= HIndex + i + 1) {
				return (uintE)(HIndex + i + 1);
			}else{
				return HIndex;
			}
		});
		uintE NHIndex = pbbs::reduce(worksOrNo,pbbs::maxm<uintE>());
		auto allH = C.A[HIndex].entries();
		par_for(0,NHIndex - (HIndex - cNs[NHIndex + 1]),[&](size_t i) {
			B.insert(allH[i]);
		});
		BSize = NHIndex - (HIndex - cNs[NHIndex + 1]);
		HIndex = NHIndex;
	}
};