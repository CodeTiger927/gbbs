#include "h-index.h"

#include <math.h>
#include <fstream>

int main() {

  // TODO: use integrated ligra main
	// e.g.,
	//template <class Graph>
  //double AppKCore_runner(Graph& GA, commandLine P) {
  //}
	//generate_symmetric_main(AppKCore_runner, false);
	graph G;


	G.insertV(0);
	G.insertV(1);
	G.insertV(2);
	G.insertV(3);
	cout << G.connect(0,1) << endl;
	cout << G.connect(0,2) << endl;
	cout << G.connect(0,3) << endl;
	cout << G.connect(1,2) << endl;
	cout << G.connect(1,3) << endl;
	cout << G.connect(2,3) << endl;





	return 0;
}