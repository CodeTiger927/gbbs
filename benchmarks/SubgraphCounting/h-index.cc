#include "h-index.h"

int main() {

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