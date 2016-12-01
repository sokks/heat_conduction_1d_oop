#include "HeatEquation.h"



HeatEquation::HeatEquation()
{
	L = 0.0;
}

HeatEquation::HeatEquation(double _L, double _lambda, double _ro, double _c, double _tLeft, double _tRight)
{
	L = _L;
	x_step = L / (X_STEPS - 1);
	lambda = _lambda;
	ro = _ro;
	c = _c;
	a_sqr = lambda / (ro * c);
	tLeft = _tLeft;
	tRight = _tRight;
	currentTime = 0;
	startTemperature = vector<double>(X_STEPS);
	for (int i = 0; i < X_STEPS; i++) {
		startTemperature[i] = fStart(x_step * i, L);
	}
	previousTemperature = vector<double>(X_STEPS, 0.0);
	currentTemperature = startTemperature;
}

HeatEquation::HeatEquation(const HeatEquation & anotherEquation)
{
	X_STEPS = anotherEquation.X_STEPS;
	TIME_STEPS = anotherEquation.TIME_STEPS;
	L = anotherEquation.L;
	lambda = anotherEquation.lambda;
	ro = anotherEquation.ro;
	c = anotherEquation.c;
	a_sqr = anotherEquation.a_sqr;
	x_step = anotherEquation.x_step;
	time_step = anotherEquation.time_step;
	tLeft = anotherEquation.tLeft;
	tRight = anotherEquation.tRight;
	startTemperature = anotherEquation.startTemperature;
	previousTemperature = anotherEquation.previousTemperature;
	currentTemperature = anotherEquation.currentTemperature;
	currentTime = anotherEquation.currentTime;
}

HeatEquation & HeatEquation::operator=(const HeatEquation & anotherEquation)
{
	X_STEPS = anotherEquation.X_STEPS;
	TIME_STEPS = anotherEquation.TIME_STEPS;
	L = anotherEquation.L;
	lambda = anotherEquation.lambda;
	ro = anotherEquation.ro;
	c = anotherEquation.c;
	a_sqr = anotherEquation.a_sqr;
	x_step = anotherEquation.x_step;
	time_step = anotherEquation.time_step;
	tLeft = anotherEquation.tLeft;
	tRight = anotherEquation.tRight;
	startTemperature = anotherEquation.startTemperature;
	previousTemperature = anotherEquation.previousTemperature;
	currentTemperature = anotherEquation.currentTemperature;
	currentTime = anotherEquation.currentTime;

	return *this;
}


HeatEquation::~HeatEquation()
{
}

void HeatEquation::SetXSteps(int newSteps)
{
	X_STEPS = newSteps;
	x_step = L / (X_STEPS - 1);
	startTemperature.resize(X_STEPS);
	for (int i = 0; i < X_STEPS; i++) {
		startTemperature[i] = fStart(x_step * i, L);
	}
	currentTime = 0.0;
	previousTemperature.clear();
	currentTemperature = startTemperature;
}

void HeatEquation::SetTimeSteps(int newTimeSteps)
{
	TIME_STEPS = newTimeSteps;
}

void HeatEquation::SetTimeStep(double newVal)
{
	time_step = newVal;
	currentTime = 0.0;
	previousTemperature.clear();
	currentTemperature = startTemperature;
}

vector<double> HeatEquation::getPreviousTemperature()
{
	return previousTemperature;
}

vector<double> HeatEquation::getCurrentTemperature()
{
	return currentTemperature;
}

double HeatEquation::getTime()
{
	return currentTime;
}

void HeatEquation::solve(double tEnd)
{
	if (tEnd - time_step * TIME_STEPS > EPS) {
		SetTimeStep(tEnd / TIME_STEPS);
	}
	for (int i = 0; i < TIME_STEPS; i++) {
		doTimeStep();
	}
}

void HeatEquation::doTimeStep()
{
	currentTime += time_step;
	previousTemperature = currentTemperature;
	double A, B, C;
	double tmp = ro * c / time_step;
	A = C = lambda / (x_step * x_step);
	B = 2 * lambda / (x_step * x_step) + tmp;
	MatrixDiag S(X_STEPS, -B, A, C);
	vector<double> F(X_STEPS);
	for (int i = 0; i < X_STEPS; i++) {
		F[i] = -tmp * previousTemperature[i];
	}
	currentTemperature = S.sweep(F);
}

double fStart(double x, double h)
{
	return 2 * sin(M_PI * x / h);
}

double testRes(double x, double t, double a, double h)
{
	return 2 * sin(M_PI * x / h) * exp(-a * M_PI * M_PI / (h * h) * t);
}

TestHeatEquation::TestHeatEquation(double _L, double _lambda, double _ro, double _c, double _tLeft, double _tRight):
									HeatEquation(_L, _lambda, _ro, _c, _tLeft, _tRight)
{
	previousTemperatureAnalytic = vector<double>(X_STEPS, 0.0);
	currentTemperatureAnalytic = startTemperature;
}

TestHeatEquation::TestHeatEquation(const TestHeatEquation & anotherEquation)
{
	X_STEPS = anotherEquation.X_STEPS;
	TIME_STEPS = anotherEquation.TIME_STEPS;
	L = anotherEquation.L;
	lambda = anotherEquation.lambda;
	ro = anotherEquation.ro;
	c = anotherEquation.c;
	a_sqr = anotherEquation.a_sqr;
	x_step = anotherEquation.x_step;
	time_step = anotherEquation.time_step;
	tLeft = anotherEquation.tLeft;
	tRight = anotherEquation.tRight;
	startTemperature = anotherEquation.startTemperature;
	previousTemperature = anotherEquation.previousTemperature;
	currentTemperature = anotherEquation.currentTemperature;
	previousTemperatureAnalytic = anotherEquation.previousTemperatureAnalytic;
	currentTemperatureAnalytic = anotherEquation.currentTemperatureAnalytic;
	currentTime = anotherEquation.currentTime;
}

TestHeatEquation & TestHeatEquation::operator=(const TestHeatEquation & anotherEquation)
{
	X_STEPS = anotherEquation.X_STEPS;
	TIME_STEPS = anotherEquation.TIME_STEPS;
	L = anotherEquation.L;
	lambda = anotherEquation.lambda;
	ro = anotherEquation.ro;
	c = anotherEquation.c;
	a_sqr = anotherEquation.a_sqr;
	x_step = anotherEquation.x_step;
	time_step = anotherEquation.time_step;
	tLeft = anotherEquation.tLeft;
	tRight = anotherEquation.tRight;
	startTemperature = anotherEquation.startTemperature;
	previousTemperature = anotherEquation.previousTemperature;
	currentTemperature = anotherEquation.currentTemperature;
	previousTemperatureAnalytic = anotherEquation.previousTemperatureAnalytic;
	currentTemperatureAnalytic = anotherEquation.currentTemperatureAnalytic;
	currentTime = anotherEquation.currentTime;

	return *this;
}

void TestHeatEquation::presolve(double tEnd, bool check, double eps)
{
	if (tEnd - time_step * TIME_STEPS > EPS) {
		SetTimeStep(tEnd / TIME_STEPS);
	}
	for (int i = 0; i < TIME_STEPS; i++) {
		doTestStep();
		if (check) {
			compare(eps);
		}
	}
}

void TestHeatEquation::doTestStep()
{
	currentTime += time_step;
	previousTemperatureAnalytic = currentTemperatureAnalytic;
	for (int i = 0; i < X_STEPS; i++) {
		currentTemperatureAnalytic[i] = testRes(x_step * i, currentTime, a_sqr, L);
	}
}

bool TestHeatEquation::compare(double eps)
{
	for (int i = 0; i < X_STEPS; i++) {
		if (currentTemperature[i] - currentTemperatureAnalytic[i] > eps)
			return false;
	}
	return true;
}
