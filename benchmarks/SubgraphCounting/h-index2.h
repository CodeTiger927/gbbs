// modification
#include <iostream>
#include <map>
#include <set>
#include <limits>
#include "ligra/ligra.h"
#include "ligra/pbbslib/sparse_table.h"

// For hashing sparse table
struct hash_uintE {
	inline size_t operator () (const uintE & a) {return pbbs::hash64_2(a);}
};

struct HSet {
	// Empty value for uintE is UINT_E_MAX
	// All vertices
	sparse_table<uintE, pbbs::empty, hash_uintE> S = 
		make_sparse_table<uintE, pbbs::empty, hash_uintE>
		(1000, std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());
	// Maps each vertex to its degree
	sparse_table<uintE, uintE, hash_uintE> f =
		make_sparse_table(1000, std::make_tuple(UINT_E_MAX, UINT_E_MAX), hash_uintE());
	// ADDDDDD COMMMMENNTSS!!!!!!!!!!!!!!!!!!!DJKDS:LKJFDS:LKJSDF:KJDSF:KJSDF
	sparse_table<uintE, pbbs::empty, hash_uintE> H = 
		make_sparse_table<uintE, pbbs::empty, hash_uintE>
		(1000, std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());
	// ADDDDDD COMMMMENNTSS!!!!!!!!!!!!!!!!!!!!!DJKDS:LKJFDS:LKJSDF:KJDSF:KJSDF
	sparse_table<uintE, pbbs::empty, hash_uintE> B = 
		make_sparse_table<uintE, pbbs::empty, hash_uintE>
		(1000, std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE());
	// ADDDDDD COMMMMENNTSS!!!!!!!!!!!!!!!!!!!!!!!DJKDS:LKJFDS:LKJSDF:KJDSF:KJSDF
	sparse_table<uintE, sparse_table<uintE, pbbs::empty, hash_uintE>, hash_uintE> C = 
		make_sparse_table<uintE, sparse_table<uintE, pbbs::empty, hash_uintE>, hash_uintE>
		(1000, std::make_tuple(UINT_E_MAX, make_sparse_table<uintE, pbbs::empty, hash_uintE>
			(0, std::make_tuple(UINT_E_MAX, pbbs::empty()), hash_uintE())), hash_uintE());



	uintE insert(uintE v, uintE deg) {
		//e = v and x = deg
		f.insert(std::make_tuple(v, deg));
		C.find(deg, make_sparse_table<uintE, pbbs::empty, hash_uintE>).insert(v); // FIX!!!!!!!DFJKHSJ:DKJ:HSDFKJHSDF::KJDSF:JSDFKJ

		if(deg > H.size()) {
/*
			H.insert(e);
			if(B.size() == 0) {
				if(C[H.size() - 1].size() == 0) {
					B = C[H.size()];
					C[H.size()] = std::set<uintE>();
				}
			}else{
				uintE cur = *(B.begin());
				B.erase(B.begin());
				H.erase(cur);
				C[H.size()].insert(cur);
			}
*/
		}

		//return H.size();

		return 0;
	}
	
	uintE erase(uintE e) {
/*
		if(B.find(e) != B.end()) {
			B.erase(e);
		}else{
			C[f[e]].erase(e);
		}
		if(H.find(e) != H.end()) {
			H.erase(e);
			if(C[H.size() + 1].size() == 0) {
				C[H.size() + 1] = B;
				B = std::set<uintE>();
			}else{
				uintE cur = *(C[H.size() + 1].begin());
				C[H.size() + 1].erase(cur);
				B.insert(cur);
				H.insert(cur);
			}
		}

		return H.size();
*/
		return 0;
	}
	
	uintE change(uintE e, uintE x) {
		erase(e);
		return insert(e,x);
	}
};

