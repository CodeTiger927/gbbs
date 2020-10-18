using namespace std;

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <assert.h>

vector<long long> edges[500000];

int main() {
	ifstream in("../graph_test_small.txt");
	ofstream out("../graph_test_small2.txt");
	for(int i = 0;i < 1049866;++i) {
		long long a,b;
		in >> a >> b;
		edges[a].push_back(b);
		edges[b].push_back(a);
	}
	out << 425876 << " " << 1049866 << endl;
	long long sum = 0;
	for(int i = 0;i < 500000;++i) {
		out << sum << "\n";
		sum += edges[i].size();
	}
	for(int i = 0;i < 500000;++i) {
		for(int j = 0;j < edges[i].size();++j) {
			out << edges[i][j] << "\n";
		}
	}
	in.close();
	out.close();
	return 0;
}
