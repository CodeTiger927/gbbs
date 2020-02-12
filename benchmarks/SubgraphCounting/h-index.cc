using namespace std;

#include "h-index.h"

int main() {
	HSet h;
	cout << h.insert(1,4) << endl;
	cout << h.insert(2,3) << endl;
	cout << h.insert(3,4) << endl;
	cout << h.insert(4,4) << endl;
	cout << h.erase(2) << endl;
	cout << h.insert(2,4) << endl;
	cout << h.erase(2) << endl;
	cout << h.insert(2,4) << endl;
	cout << h.change(2,3) << endl;
	cout << h.change(2,4) << endl;
	return 0;
}