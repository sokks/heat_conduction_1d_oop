#include "HeatEquation.h"
#include <iostream>

using std::cout;
using std::endl;

int main() {
	//  сталь
	/*double lambda = 46.0;
	double ro = 7800.0;
	double c = 460.0;
	int t_end = 120;
	double L = 0.1;
	char output[] = "output.txt";*/
	HeatEquation Eq1(1, 46.0, 7800.0, 460.0, 0.0, 0.0);
	Eq1.solve(20.0);
	vector<double> res = Eq1.getCurrentTemperature();
	for (vector<double>::iterator it = res.begin(); it != res.end(); ++it) {
		cout << *it << " ";
	}
	cout << endl;
	return 0;
}