#include "MatrixDiag.h"



MatrixDiag::MatrixDiag()
{
	size = 0;
}

MatrixDiag::MatrixDiag(int n, double diagElement, double upElement, double downElement)
{
	size = n;
	mainDiag = vector<double>(n, diagElement);
	upDiag = vector<double>(n, upElement);
	downDiag = vector<double>(n, downElement);
}

MatrixDiag::MatrixDiag(const MatrixDiag & anotherMatrixDiag)
{
	size = anotherMatrixDiag.size;
	mainDiag = anotherMatrixDiag.mainDiag;
	upDiag = anotherMatrixDiag.upDiag;
	downDiag = anotherMatrixDiag.downDiag;
}

MatrixDiag & MatrixDiag::operator=(const MatrixDiag & anotherMatrixDiag)
{
	size = anotherMatrixDiag.size;
	mainDiag = anotherMatrixDiag.mainDiag;
	upDiag = anotherMatrixDiag.upDiag;
	downDiag = anotherMatrixDiag.downDiag;
	return *this;
}


MatrixDiag::~MatrixDiag()
{
	//  standart delete
	//  vector destructor works
}

bool MatrixDiag::isEmpty()
{
	return (size != 0);
}

vector<double> MatrixDiag::sweep(vector<double> F)
{
	vector<double> X(size, 0.0);
	vector<double> alfa(size, 0.0);
	vector<double> beta(size, 0.0);
	__int64 ctr1 = 0, ctr2 = 0, freq = 0;
	QueryPerformanceCounter((LARGE_INTEGER *)&ctr1);
	alfa[1] = -upDiag[0] / mainDiag[0];
	beta[1] = F[0] / mainDiag[0];
	for (int i = 2; i < size; i++) {
		alfa[i] = -upDiag[i - 1] / (downDiag[i - 1] * alfa[i - 1] + mainDiag[i - 1]);
		beta[i] = (F[i - 1] - downDiag[i - 1] * beta[i - 1]) / (downDiag[i - 1] * alfa[i - 1] + mainDiag[i - 1]);
	}

	X[size - 1] = (F[size - 1] - downDiag[size - 1] * beta[size - 1]) / (downDiag[size - 1] * alfa[size - 1] + mainDiag[size - 1]);
	for (int i = size - 2; i >= 0; i--) {
		X[i] = alfa[i + 1] * X[i + 1] + beta[i + 1];
	}
	QueryPerformanceCounter((LARGE_INTEGER *)&ctr2);
	QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
	double wtime((ctr2 - ctr1) * 1.0 / freq);
	ofstream fout("wtime1.txt", std::fstream::app);
	fout << wtime << endl;
	fout.close();
	return X;
}

vector<double> MatrixDiag::sweepOMP(vector<double> F, int nOfThreads)
{
	int curNOfThreads = omp_get_num_threads();
	omp_set_num_threads(nOfThreads);
	vector<double> X(size);

	

	omp_set_num_threads(curNOfThreads);
	return X;
}

vector<double> MatrixDiag::sweepOMP1(vector<double> F)
{
	int curNOfThreads = omp_get_num_threads();
	omp_set_num_threads(2);
	double wtime0 = 0.0, wtime1 = 0.0;
	vector<double> X(size, 0.0);
	vector<double> alfa(size, 0.0);
	vector<double> beta(size, 0.0);
	vector<double> ksi(size, 0.0);
	vector<double> eta(size, 0.0);
	int p = size / 2;

#pragma omp parallel shared(X, alfa, beta, ksi, eta, p)
{
#pragma omp sections firstprivate(wtime0, wtime1)
	{
#pragma omp section
		{
			//  thread 1
			wtime0 = omp_get_wtime();
			alfa[1] = -upDiag[0] / mainDiag[0];
			beta[1] = F[0] / mainDiag[0];
			for (int i = 2; i <= p + 1; i++) {
				alfa[i] = -upDiag[i - 1] / (downDiag[i - 1] * alfa[i - 1] + mainDiag[i - 1]);
				beta[i] = (F[i - 1] - downDiag[i - 1] * beta[i - 1]) / (downDiag[i - 1] * alfa[i - 1] + mainDiag[i - 1]);
			}
			wtime0 = omp_get_wtime() - wtime0;
		}
#pragma omp section
		{
			//  thread 2
			wtime0 = omp_get_wtime();
			ksi[size - 1] = -downDiag[size - 1] / mainDiag[size - 1];
			eta[size - 1] = F[size - 1] / mainDiag[size - 1];
			for (int i = size - 2; i > p; i--) {
				ksi[i] = -downDiag[i] / (upDiag[i] * ksi[i + 1] + mainDiag[i]);
				eta[i] = (F[i] - upDiag[i] * eta[i + 1]) / (upDiag[i] * ksi[i + 1] + mainDiag[i]);
			}
			wtime0 = omp_get_wtime() - wtime0;
		}
	}
	//  one thread
#pragma omp barrier
#pragma omp single 
	{
		X[p] = (beta[p + 1] + alfa[p + 1] * eta[p + 1]) / (1 - alfa[p + 1] * ksi[p + 1]);
	}
#pragma omp barrier
#pragma omp sections
	{
#pragma omp section
		{
			//  thread 1
			wtime1 = omp_get_wtime();
			for (int i = p - 1; i >= 0; i--) {
				X[i] = alfa[i + 1] * X[i + 1] + beta[i + 1];
			}
			wtime1 = omp_get_wtime() - wtime1;
		}
#pragma omp section
		{
			//  thread 2
			wtime1 = omp_get_wtime();
			for (int i = p + 1; i < size; i++) {
				X[i] = ksi[i] * X[i - 1] + eta[i];
			}
			wtime1 = omp_get_wtime() - wtime1;
		}
	}
}
	omp_set_num_threads(curNOfThreads);
	ofstream fout("wtime2.txt", std::fstream::app);
	fout << wtime0 + wtime1 << endl;
	fout.close();
	return X;
}


