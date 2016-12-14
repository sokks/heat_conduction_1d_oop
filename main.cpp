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
	char out1[] = "output0.txt";
	char out2[] = "outputtest0.txt";
	TestHeatEquation Eq1(0.1, 46.0, 7800.0, 460.0, 0.0, 0.0);
	//TestHeatEquation Eq1(1.0, 1.0, 1.0, 1.0, 0.0, 0.0);
	cout << Eq1.getCurrentTemperature() << endl << endl;
	cout << Eq1.getPreviousTemperature() << endl << endl;
	//double wtime2 = Eq1.solveOMP1(out1, 1.0);
	double wtime1 = Eq1.solve(out1, 1.0);
	cout << "Time of solving: " << wtime1 << endl;
	//cout << "Time of OMP solving: " << wtime2 << endl;
	//ofstream fout("testres3.txt", std::fstream::app);
	//fout << wtime1 / wtime2 << endl;
	//fout.close();
	Eq1.presolve(out2, 1.0);
	cout << endl;
	
	//Eq1.doOMPTimeStep();
	return 0;
}