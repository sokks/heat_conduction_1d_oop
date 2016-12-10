#include "HeatEquation.h"
#include <iostream>

using std::cout;
using std::endl;

int main() {
	//omp_set_dynamic(1);
	//cout << omp_get_dynamic() << endl;

	//  �����
	/*double lambda = 46.0;
	double ro = 7800.0;
	double c = 460.0;
	int t_end = 120;
	double L = 0.1;
	char output[] = "output.txt";*/
	std::string out1("output.txt");
	std::string out2("outputtest.txt");
	TestHeatEquation Eq1(0.1, 46.0, 7800.0, 460.0, 0.0, 0.0);
	//TestHeatEquation Eq1(1.0, 1.0, 1.0, 1.0, 0.0, 0.0);
	cout << "Time of solving:" << Eq1.solve(out1, 1.0) << endl;
	cout << "Time of OMP solving:" << Eq1.solveOMP(out1, 1.0) << endl;
	Eq1.presolve(out2, 1.0);
	cout << endl;
	return 0;
}