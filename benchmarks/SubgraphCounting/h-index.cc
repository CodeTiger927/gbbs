#include "h-index.h"

int main() {

	// graph G;


	// G.insertV(0);
	// G.insertV(1);
	// G.insertV(2);
	// G.insertV(3);
	// cout << G.connect(0,1) << endl;
	// cout << G.connect(0,2) << endl;
	// cout << G.connect(0,3) << endl;
	// cout << G.connect(1,2) << endl;
	// cout << G.connect(1,3) << endl;
	// cout << G.connect(2,3) << endl;

	dynamic_vertex_data dvd;
	dvd.degree = 3;
  dvd.neighbors = make_sparse_table<uintE,bool,hash_uintE>(1000,std::make_tuple(UINT_E_MAX,false),hash_uintE());
  dvd.neighbors.insert(std::make_tuple(3,true));
  dvd.neighbors.insert(std::make_tuple(5,true));
  dvd.neighbors.insert(std::make_tuple(7,true));




	dynamic_symmetric_vertex<uintE> sv = dynamic_symmetric_vertex<uintE>(dvd);

  cout << sv.getInNeighbor(2) << endl;



	return 0;
}