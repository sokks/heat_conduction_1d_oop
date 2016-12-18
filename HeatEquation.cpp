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
	startTemperature = (double *)malloc(sizeof(double) * X_STEPS);
	for (int i = 0; i < X_STEPS; i++) {
		startTemperature[i] = fStart(x_step * i, L);  //f1(x_step * i)
	}
	previousTemperature = (double *)malloc(sizeof(double) * X_STEPS);
	for (int i = 0; i < X_STEPS; i++)
		previousTemperature[i] = 0.0;
	currentTemperature = (double *)malloc(sizeof(double) * X_STEPS);
	for (int i = 0; i < X_STEPS; i++)
		currentTemperature[i] = startTemperature[i];
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
	startTemperature = (double *)malloc(sizeof(double) * X_STEPS); 
	for (int i = 0; i < X_STEPS; i++)
		startTemperature[i] = anotherEquation.startTemperature[i];
	previousTemperature = (double *)malloc(sizeof(double) * X_STEPS);
	for (int i = 0; i < X_STEPS; i++)
		previousTemperature[i] = anotherEquation.previousTemperature[i];
	currentTemperature = (double *)malloc(sizeof(double) * X_STEPS);
	for (int i = 0; i < X_STEPS; i++)
		currentTemperature[i] = anotherEquation.currentTemperature[i];
	currentTime = anotherEquation.currentTime;
}

HeatEquation & HeatEquation::operator=(const HeatEquation & anotherEquation)
{
	if (this == &anotherEquation)
		return *this;

	if (anotherEquation.X_STEPS != X_STEPS) { //?
		X_STEPS = anotherEquation.X_STEPS;
		free(startTemperature);
		startTemperature = (double *)malloc(sizeof(double) * X_STEPS);
		free(previousTemperature);
		previousTemperature = (double *)malloc(sizeof(double) * X_STEPS);
		free(currentTemperature);
		currentTemperature = (double *)malloc(sizeof(double) * X_STEPS);
	}

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
	for (int i = 0; i < X_STEPS; i++)
		startTemperature[i] = anotherEquation.startTemperature[i];
	for (int i = 0; i < X_STEPS; i++)
		previousTemperature[i] = anotherEquation.previousTemperature[i];
	for (int i = 0; i < X_STEPS; i++)
		currentTemperature[i] = anotherEquation.currentTemperature[i];
	currentTime = anotherEquation.currentTime;

	return *this;
}


HeatEquation::~HeatEquation()
{
	if (startTemperature != 0) {
		free(startTemperature);
		free(previousTemperature);
		free(currentTemperature);
	}
}

void HeatEquation::SetXSteps(int newSteps)
{
	X_STEPS = newSteps;
	x_step = L / (X_STEPS - 1);
	free(startTemperature);
	startTemperature = (double *)malloc(sizeof(double) * X_STEPS);
	for (int i = 0; i < X_STEPS; i++) {
		startTemperature[i] = fStart(x_step * i, L); //f1(x_step * i)
	}
	currentTime = 0.0;
	free(previousTemperature);
	previousTemperature = (double *)malloc(sizeof(double) * X_STEPS);
	for (int i = 0; i < X_STEPS; i++)
		previousTemperature[i] = 0.0;
	currentTemperature = (double *)malloc(sizeof(double) * X_STEPS);
	for (int i = 0; i < X_STEPS; i++)
		currentTemperature[i] = startTemperature[i];
}

void HeatEquation::SetTimeSteps(int newTimeSteps)
{
	TIME_STEPS = newTimeSteps;
}

void HeatEquation::SetTimeStep(double newVal)
{
	time_step = newVal;
	currentTime = 0.0;
	free(previousTemperature);
	for (int i = 0; i < X_STEPS; i++) {
		previousTemperature[i] = 0.0;
		currentTemperature[i] = startTemperature[i];
	}
}

double * HeatEquation::getPreviousTemperature()
{
	return previousTemperature;
}

double * HeatEquation::getCurrentTemperature()
{
	return currentTemperature;
}

double HeatEquation::getTime()
{
	return currentTime;
}

double HeatEquation::solve(char *filename, double tEnd)
{
	if (tEnd - time_step * TIME_STEPS > EPS) {
		SetTimeStep(tEnd / TIME_STEPS);
	}
	double runtime = 0.0;
	currentTime = 0;
	for (int i = 0; i < X_STEPS; i++)
		currentTemperature[i] = startTemperature[i];
	ofstream fout;
	fout.open(filename);
	fout << TIME_STEPS << std::endl;
	for (int i = 0; i < TIME_STEPS; i++) {
		for (int j = 0; j < X_STEPS; j++) {
			fout << currentTemperature[j] << " ";
		}
		fout << std::endl;
		runtime += doTimeStep();
	} // не вывожу последнюю температуру?
	fout.close();
	return runtime;
}

double HeatEquation::doTimeStep()
{
	currentTime += time_step;
	__int64 ctr1 = 0, ctr2 = 0, freq = 0;
	QueryPerformanceCounter((LARGE_INTEGER *)&ctr1);
	for (int i = 0; i < X_STEPS; i++)
		previousTemperature[i] = currentTemperature[i];
	free(currentTemperature);
	currentTemperature = 0;
	double A, B, C;
	double tmp = ro * c / time_step;
	A = C = lambda / (x_step * x_step);
	B = 2 * lambda / (x_step * x_step) + tmp;
	MatrixDiag S(X_STEPS, -B, A, C);
	double *F = (double *)malloc(sizeof(double) * X_STEPS);
	for (int i = 0; i < X_STEPS; i++) {
		F[i] = -tmp * previousTemperature[i];
	}
	currentTemperature = S.sweep(F);
	free(F);
	QueryPerformanceCounter((LARGE_INTEGER *)&ctr2);
	QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
	return ((ctr2 - ctr1) * 1.0 / freq); //  время в микросекундах
}

double HeatEquation::doOMPTimeStep()
{
	currentTime += time_step;
	__int64 ctr1 = 0, ctr2 = 0, freq = 0;
	QueryPerformanceCounter((LARGE_INTEGER *)&ctr1);

	for (int i = 0; i < X_STEPS; i++)
		previousTemperature[i] = currentTemperature[i];
	free(currentTemperature);
	currentTemperature = 0;

	double A, B, C;
	double tmp;
	tmp = ro * c / time_step;
	A = C = lambda / (x_step * x_step);
	B = 2 * lambda / (x_step * x_step) + tmp;
	MatrixDiag S(X_STEPS, -B, A, C, false, true);
	double *F = (double *)malloc(sizeof(double) * X_STEPS);
#pragma omp parallel for schedule (guided, 100)
	for (int i = 0; i < X_STEPS; i++) {
		F[i] = -tmp * previousTemperature[i];
	}
	currentTemperature = S.sweepOMP(F, 4);
	free(F);
	QueryPerformanceCounter((LARGE_INTEGER *)&ctr2);
	QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
	return ((ctr2 - ctr1) * 1.0 / freq); //  время в микросекундах
}

double HeatEquation::doOMPTimeStep1()
{
	currentTime += time_step;
	__int64 ctr1 = 0, ctr2 = 0, freq = 0;
	QueryPerformanceCounter((LARGE_INTEGER *)&ctr1);
	for (int i = 0; i < X_STEPS; i++)
		previousTemperature[i] = currentTemperature[i];
	free(currentTemperature);
	currentTemperature = 0;
	double A, B, C;
	double tmp;
	tmp = ro * c / time_step;
	A = C = lambda / (x_step * x_step);
	B = 2 * lambda / (x_step * x_step) + tmp;
	MatrixDiag S(X_STEPS, -B, A, C);
	double *F = (double *)malloc(sizeof(double) * X_STEPS);

#pragma omp parallel for schedule (static)
	for (int i = 0; i < X_STEPS; i++) {
		F[i] = -tmp * previousTemperature[i];
	}
	
	currentTemperature = S.sweepOMP1(F); // BAD! maybe safe pointer class?
	free(F);
	QueryPerformanceCounter((LARGE_INTEGER *)&ctr2);
	QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
	return ((ctr2 - ctr1) * 1.0 / freq); //  время в микросекундах
}

double HeatEquation::solveOMP1(char *filename, double tEnd)
{
	if (tEnd - time_step * TIME_STEPS > EPS) {
		SetTimeStep(tEnd / TIME_STEPS);
	}
	double runtime = 0.0;
	currentTime = 0;
	for (int i = 0; i < X_STEPS; i++)
		currentTemperature[i] = startTemperature[i];
	ofstream fout;
	fout.open(filename);
	fout << TIME_STEPS << std::endl;
	for (int i = 0; i < TIME_STEPS; i++) {
		for (int j = 0; j < X_STEPS; j++) {
			fout << currentTemperature[j] << " ";
		}
		fout << std::endl;
		runtime += doOMPTimeStep1();
	}
	fout.close();
	return runtime;
}

double HeatEquation::solveOMP(char *filename, double tEnd)
{
	if (tEnd - time_step * TIME_STEPS > EPS) {
		SetTimeStep(tEnd / TIME_STEPS);
	}
	double runtime = 0.0;
	currentTime = 0;
	for (int i = 0; i < X_STEPS; i++)
		currentTemperature[i] = startTemperature[i];
	ofstream fout;
	fout.open(filename);
	fout << TIME_STEPS << std::endl;
	for (int i = 0; i < TIME_STEPS; i++) {
		for (int j = 0; j < X_STEPS; j++) {
			fout << currentTemperature[j] << " ";
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
	previousTemperatureAnalytic = (double *)malloc(sizeof(double) * X_STEPS);
	for (int i = 0; i < X_STEPS; i++)
		previousTemperatureAnalytic[i] = 0.0;
	currentTemperatureAnalytic = (double *)malloc(sizeof(double) * X_STEPS);
	for (int i = 0; i < X_STEPS; i++)
		currentTemperatureAnalytic[i] = startTemperature[i];
}

TestHeatEquation::TestHeatEquation(const TestHeatEquation & anotherEquation) : 
									HeatEquation((HeatEquation)anotherEquation)
{
	/*
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
	*/
	previousTemperatureAnalytic = (double *)malloc(sizeof(double) * X_STEPS);
	for (int i = 0; i < X_STEPS; i++)
		previousTemperatureAnalytic[i] = anotherEquation.previousTemperatureAnalytic[i];
	currentTemperatureAnalytic = (double *)malloc(sizeof(double) * X_STEPS);
	for (int i = 0; i < X_STEPS; i++)
		currentTemperatureAnalytic[i] = anotherEquation.currentTemperatureAnalytic[i];
}

TestHeatEquation & TestHeatEquation::operator=(const TestHeatEquation & anotherEquation)
{
	if (&anotherEquation == this) {
		return *this;
	}

	(HeatEquation)(*this) = (HeatEquation)anotherEquation;
	
	if (X_STEPS != anotherEquation.X_STEPS) {
		free(previousTemperatureAnalytic);
		previousTemperatureAnalytic = (double *)malloc(sizeof(double) * X_STEPS); 
		for (int i = 0; i < X_STEPS; i++) 
			previousTemperatureAnalytic[i] = anotherEquation.previousTemperatureAnalytic[i];
		free(currentTemperatureAnalytic);
		currentTemperatureAnalytic = (double *)malloc(sizeof(double) * X_STEPS);
		for (int i = 0; i < X_STEPS; i++)
			currentTemperatureAnalytic[i] = anotherEquation.currentTemperatureAnalytic[i];
	}
	return *this;
}

void TestHeatEquation::presolve(char *filename, double tEnd, bool check, double eps)
{
	if (tEnd - time_step * TIME_STEPS > EPS) {
		SetTimeStep(tEnd / TIME_STEPS);
	}
	currentTime = 0;
	for (int i = 0; i < X_STEPS; i++)
		currentTemperatureAnalytic[i] = startTemperature[i];
	ofstream fout;
	fout.open(filename);
	fout << TIME_STEPS << std::endl;
	for (int i = 0; i < TIME_STEPS; i++) {
		for (int j = 0; j < X_STEPS; j++) {
			fout << currentTemperatureAnalytic[j] << " ";
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
	for (int i = 0; i < X_STEPS; i++) {
		previousTemperatureAnalytic[i] = currentTemperatureAnalytic[i];
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
