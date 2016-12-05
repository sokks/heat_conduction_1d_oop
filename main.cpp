#include "HeatEquation.h"
#include <iostream>

using std::cout;
using std::endl;

int main() {
	//omp_set_dynamic(1);
	//cout << omp_get_dynamic() << endl;

	//  сталь
	/*double lambda = 46.0;
	double ro = 7800.0;
	double c = 460.0;
	int t_end = 120;
	double L = 0.1;
	char output[] = "output.txt";*/
	std::string out1("output.txt");
	std::string out2("outputtest.txt");
	TestHeatEquation Eq1(0.1, 46.0, 7800.0, 460.0, 0.0, 0.0);
	Eq1.solve(out1, 20.0);
	Eq1.presolve(out2, 20.0);
	cout << endl;
	return 0;
}