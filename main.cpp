#include "HeatEquation.h"
#include <iostream>

using std::cout;
using std::cin;
using std::endl;

int main() {
 
	/* сталь
	double lambda = 46.0;
	double ro = 7800.0;
	double c = 460.0;
	int t_end = 120;
	double L = 0.1;
	*/

	std::string out1("output_i.txt");
	std::string out2("output_e.txt");
	std::string out3("output_t.txt");

	double L = 1.0;
	double a = 1.0;
	int N = 100;
	int mod = 1;

	cout << "L = ";
	cin >> L;
	cout << "a = ";
	cin >> a;
	cout << "merge_mod: ";
	cin >> mod;

	TestHeatEquation Eq1(L, a, mod);

	cout << "X_STEPS = ";
	cin >> N;
	Eq1.SetXSteps(N);
	//vector<double> cur(N);

	Eq1.solve_implicit(out1, 0.5);
	Eq1.solve_explicit(out2, 0.5);

	//Eq1.presolve(out3, 0.5);

	cout << "done" << endl;
	
	return 0;
}