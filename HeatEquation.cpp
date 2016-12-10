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
		startTemperature[i] = fStart(x_step * i, L);  //f1(x_step * i)
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
		startTemperature[i] = fStart(x_step * i, L); //f1(x_step * i)
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

double HeatEquation::solve(std::string filename, double tEnd)
{
	if (tEnd - time_step * TIME_STEPS > EPS) {
		SetTimeStep(tEnd / TIME_STEPS);
	}
	double runtime = 0.0;
	currentTime = 0;
	currentTemperature = startTemperature;
	ofstream fout;
	fout.open(filename);
	fout << TIME_STEPS << std::endl;
	for (int i = 0; i < TIME_STEPS; i++) {
		for (vector<double>::iterator it = currentTemperature.begin(); it != currentTemperature.end(); ++it) {
			fout << *it << " ";
		}
		fout << std::endl;
		runtime += doTimeStep();
	}
	fout.close();
	return runtime;
}

double HeatEquation::doTimeStep()
{
	currentTime += time_step;
	__int64 ctr1 = 0, ctr2 = 0, freq = 0;
	QueryPerformanceCounter((LARGE_INTEGER *)&ctr1);
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
	QueryPerformanceCounter((LARGE_INTEGER *)&ctr2);
	QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
	return ((ctr2 - ctr1) * 1.0 / freq); //  время в микросекундах
}

double HeatEquation::doOMPTimeStep()
{
	currentTime += time_step;
	__int64 ctr1 = 0, ctr2 = 0, freq = 0;
	QueryPerformanceCounter((LARGE_INTEGER *)&ctr1);
	previousTemperature = currentTemperature;
	double A, B, C;
	double tmp;
	tmp = ro * c / time_step;
	A = C = lambda / (x_step * x_step);
	B = 2 * lambda / (x_step * x_step) + tmp;
	MatrixDiag S(X_STEPS, -B, A, C);
	vector<double> F(X_STEPS);
#pragma omp parallel for schedule (guided, 10)
	for (int i = 0; i < X_STEPS; i++) {
		F[i] = -tmp * previousTemperature[i];
	}
	currentTemperature = S.sweepOMP1(F);
	QueryPerformanceCounter((LARGE_INTEGER *)&ctr2);
	QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
	return ((ctr2 - ctr1) * 1.0 / freq); //  время в микросекундах
}

double HeatEquation::solveOMP(std::string filename, double tEnd)
{
	if (tEnd - time_step * TIME_STEPS > EPS) {
		SetTimeStep(tEnd / TIME_STEPS);
	}
	double runtime = 0.0;
	currentTime = 0;
	currentTemperature = startTemperature;
	ofstream fout;
	fout.open(filename);
	fout << TIME_STEPS << std::endl;
	for (int i = 0; i < TIME_STEPS; i++) {
		for (vector<double>::iterator it = currentTemperature.begin(); it != currentTemperature.end(); ++it) {
			fout << *it << " ";
		}
		fout << std::endl;
		runtime += doOMPTimeStep();
	}
	fout.close();
	return runtime;
}

double fStart(double x, double h)
{
	return 2 * sin(M_PI * x / h);
}

double testRes(double x, double t, double a, double h)
{
	return 2 * sin(M_PI * x / h) * exp( (-a) * M_PI * M_PI * t / (h * h) );
}

double f1(double x)
{
	return sin(M_PI * x);
}

double u1(double x, double t)
{
	return exp(-M_PI * M_PI * t) * sin(M_PI * x);
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

void TestHeatEquation::presolve(std::string filename, double tEnd, bool check, double eps)
{
	if (tEnd - time_step * TIME_STEPS > EPS) {
		SetTimeStep(tEnd / TIME_STEPS);
	}
	currentTime = 0;
	currentTemperatureAnalytic = startTemperature;
	ofstream fout;
	fout.open(filename);
	fout << TIME_STEPS << std::endl;
	for (int i = 0; i < TIME_STEPS; i++) {
		for (vector<double>::iterator it = currentTemperatureAnalytic.begin(); it != currentTemperatureAnalytic.end(); ++it) {
			fout << *it << " ";
		}
		fout << std::endl;
		doTestStep();
		if (check) {
			compare(eps); // write somewhere
		}
	}
	fout.close();
}

void TestHeatEquation::doTestStep()
{
	currentTime += time_step;
	previousTemperatureAnalytic = currentTemperatureAnalytic;
	for (int i = 0; i < X_STEPS; i++) {
		currentTemperatureAnalytic[i] = testRes(x_step * i, currentTime, a_sqr, L);  //u1(x_step * i, currentTime)
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

double Functor1::operator()()
{
	return 0.0;
}

double Functor1::operator()(double x, double h)
{
	return 2 * sin(M_PI * x / h);
}

double Functor1::operator()(double x, double t, double a, double h)
{
	return 2 * sin(M_PI * x / h) * exp((-a) * M_PI * M_PI * t / (h * h));
}

double Functor2::operator()(double x)
{
	return sin(M_PI * x);
}

double Functor2::operator()(double x, double t)
{
	return exp(-M_PI * M_PI * t) * sin(M_PI * x);
}
