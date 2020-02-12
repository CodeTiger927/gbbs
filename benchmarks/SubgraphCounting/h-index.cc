using namespace std;

#include <iostream>
#include <map>
#include <set>

set<int> S;
map<int,int> f;
set<int> H;
map<int,int> F;
set<int> B;
map<int,set<int>> C;

int insert(int e,int x) {
	F[x] = e;
	f[e] = x;
	C[x].insert(e);
	if(x > H.size()) {
		H.insert(e);
		if(B.size() == 0) {
			if(C[H.size() - 1].size() == 0) {
				B = C[H.size()];
				C[H.size()] = set<int>();
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
	F.erase(e);
	if(B.find(e) != B.end()) {
		B.erase(e);
	}else{
		C[f[e]].erase(e);
	}
	if(H.find(e) != H.end()) {
		H.erase(e);
		if(C[H.size() + 1].size() == 0) {
			C[H.size() + 1] = B;
			B = set<int>();
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

int main() {
	cout << insert(1,4) << endl;
	cout << insert(2,3) << endl;
	cout << insert(3,4) << endl;
	cout << insert(4,4) << endl;
	cout << erase(2) << endl;
	cout << insert(2,4) << endl;
	cout << erase(2) << endl;
	cout << insert(2,4) << endl;
	cout << change(2,3) << endl;
	cout << change(2,4) << endl;
	return 0;
}