#pragma once
#include <vector>
#include <string>
#include <fstream>
#include "MatrixDiag.h"
#define M_PI 3.14159265358979323846
#define EPS 0.00001
using namespace std;

class Functor1 {
	double operator() ();
	double operator() (double x, double h);
	double operator() (double x, double t, double a, double h);
};

class Functor2 {
	double operator() (double x);
	double operator() (double x, double t);
};

//  одномерное линейное однородное уравнение теплопроводности
class HeatEquation
{
protected:
	int X_STEPS = 100, TIME_STEPS = 100;
	double L = 1.0, lambda, ro, c, a_sqr;
	double x_step, time_step = 0.1;
	double tLeft = 0.0, tRight = 0.0;
	vector<double> startTemperature, previousTemperature, currentTemperature;
	double currentTime;

	double fStart(double x);
	double mu1(double t);
	double mu2(double t);

public:
	HeatEquation();
	HeatEquation(double _L, double _lambda, double _ro, double _c, double _tLeft, double _tRight);
	HeatEquation(double _L, double _a);
	HeatEquation(const HeatEquation & anotherEquation);
	HeatEquation & operator= (const HeatEquation & anotherEquation);
	~HeatEquation();

	void SetXSteps(int newSteps); //  сбрасывает время и результаты
	void SetTimeSteps(int newTimeSteps); 
	void SetTimeStep(double newVal);  //  тоже сбрасывает

	vector<double> getPreviousTemperature();
	vector<double> getCurrentTemperature();
	double getTime();

	double solve_implicit(std::string filename, double tEnd = 10.0);
	//double solve_explicit(std::string filename, double tEnd = 10.0);
	double doTimeStep();
};

class TestHeatEquation : public HeatEquation {
	vector<double> previousTemperatureAnalytic, currentTemperatureAnalytic;
	double testRes(double x, double t);
public:
	TestHeatEquation() {}
	TestHeatEquation(double _L, double _lambda, double _ro, double _c, double _tLeft, double _tRight);
	TestHeatEquation(double _L, double _a);
	TestHeatEquation(const TestHeatEquation & anotherEquation);
	TestHeatEquation & operator= (const TestHeatEquation & anotherEquation);
	void presolve(std::string filename, double tEnd = 10.0, bool check = false, double eps = 0.001);
	void doTestStep();
	bool compare(double eps = 0.001);
};

