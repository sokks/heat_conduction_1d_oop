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

double f1(double x);
double u1(double x, double t);

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
	int X_STEPS = 1024, TIME_STEPS = 100;
	double L, lambda, ro, c, a_sqr;
	double x_step, time_step = 0.1;
	double tLeft, tRight;
	double *startTemperature, *previousTemperature, *currentTemperature;
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

	double *getPreviousTemperature();
	double *getCurrentTemperature();
	double getTime();
	int getXSteps() { return X_STEPS; }

	double solve(char *filename, double tEnd = 10.0);
	double doTimeStep();
	double doOMPTimeStep();
	double doOMPTimeStep1();
	double solveOMP1(char *filename, double tEnd = 10.0);
	double solveOMP(char *filename, double tEnd = 10.0);

};

class TestHeatEquation : public HeatEquation {
	double *previousTemperatureAnalytic, *currentTemperatureAnalytic;
public:
	TestHeatEquation() {}
	TestHeatEquation(double _L, double _lambda, double _ro, double _c, double _tLeft, double _tRight);
	TestHeatEquation(const TestHeatEquation & anotherEquation);
	TestHeatEquation & operator= (const TestHeatEquation & anotherEquation);
	void presolve(char *filename, double tEnd = 10.0, bool check = false, double eps = 0.001);
	void doTestStep();
	bool compare(double eps = 0.001);
};

