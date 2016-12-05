#pragma once
#include <vector>
#include <string>
#include <fstream>
#include "MatrixDiag.h"
#define M_PI 3.14159265358979323846
#define EPS 0.00001
using namespace std;

double fStart(double x, double h);
double testRes(double x, double t, double a, double h);

//  одномерное линейное однородное уравнение теплопроводности
class HeatEquation
{
protected:
	int X_STEPS = 101, TIME_STEPS = 100;
	double L, lambda, ro, c, a_sqr;
	double x_step, time_step = 0.1;
	double tLeft, tRight;
	vector<double> startTemperature, previousTemperature, currentTemperature;
	double currentTime;

public:
	HeatEquation();
	HeatEquation(double _L, double _lambda, double _ro, double _c, double _tLeft, double _tRight);
	HeatEquation(const HeatEquation & anotherEquation);
	HeatEquation & operator= (const HeatEquation & anotherEquation);
	~HeatEquation();

	void SetXSteps(int newSteps); //  сбрасывает время и результаты
	void SetTimeSteps(int newTimeSteps); 
	void SetTimeStep(double newVal);  //  тоже сбрасывает

	vector<double> getPreviousTemperature();
	vector<double> getCurrentTemperature();
	double getTime();

	void solve(std::string filename, double tEnd = 10.0);
	void doTimeStep();

};

class TestHeatEquation : public HeatEquation {
	vector<double> previousTemperatureAnalytic, currentTemperatureAnalytic;
public:
	TestHeatEquation() {}
	TestHeatEquation(double _L, double _lambda, double _ro, double _c, double _tLeft, double _tRight);
	TestHeatEquation(const TestHeatEquation & anotherEquation);
	TestHeatEquation & operator= (const TestHeatEquation & anotherEquation);
	void presolve(std::string filename, double tEnd = 10.0, bool check = false, double eps = 0.001);
	void doTestStep();
	bool compare(double eps = 0.001);
};

