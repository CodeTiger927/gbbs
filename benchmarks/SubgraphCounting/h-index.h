#include <iostream>
#include <map>
#include <set>
#include <limits>
#include "ligra/ligra.h"
#include "ligra/pbbslib/sparse_table.h"

struct HSet {
	std::set<int> S;
	int f[1000];
	std::set<int> H;
	std::set<int> B;
	std::set<int> C[1000];

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

