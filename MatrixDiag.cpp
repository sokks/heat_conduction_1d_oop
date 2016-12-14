#include "MatrixDiag.h"



MatrixDiag::MatrixDiag()
{
	size = 0;
	mainElem = upElem = downElem = 0;
	mainDiag = upDiag = downDiag = Vec();
	matrix = 0;
}

MatrixDiag::MatrixDiag(int n, double diagElement, double upElement, double downElement, bool stat, bool full)
{
	size = n;
	
	if (stat) {
		is_static = true;
		mainElem = diagElement;
		upElem = upElement;
		downElem = downElement;
		mainDiag = upDiag = downDiag = Vec();
	}
	else {
		mainElem = upElem = downElem = 0;
		mainDiag = Vec(size, diagElement);
		upDiag = Vec(size, upElement);
		downDiag = Vec(size, downElement);
	}

	if (full) {
		is_full = true;
		matrix = (double **)malloc(sizeof(double *) * size);
		for (int i = 0; i < size; i++)
			matrix[i] = (double *)malloc(sizeof(double) * size);
		matrix[0][0] = diagElement;
		matrix[0][1] = upElement;
		for (int i = 2; i < size; i++)
			matrix[0][i] = 0.0;
		int tmp = size - 1, tmp2 = size - 2;
		matrix[tmp][tmp] = diagElement;
		matrix[tmp][tmp2] = downElement;
		for (int i = 0; i < tmp2; i++)
			matrix[tmp][i] = 0.0;
		for (int i = 1; i < tmp; i++) {
			matrix[i][i] = diagElement;
			matrix[i][i - 1] = downElement;
			matrix[i][i + 1] = upElement;
		}
	}
	else {
		matrix = 0;
	}
}

MatrixDiag::MatrixDiag(int n, Vec diag, Vec up, Vec down)
{
	size = n;
	mainElem = upElem = downElem = 0;
	mainDiag = diag;
	upDiag = up;
	downDiag = down;
}

MatrixDiag::MatrixDiag(const MatrixDiag & anotherMatrixDiag)
{
	size = anotherMatrixDiag.size;

	if (size == 0) {
		mainElem = upElem = downElem = 0;
		mainDiag = upDiag = downDiag = Vec();
		matrix = 0;
		return;
	}

	if (anotherMatrixDiag.is_static) {
		is_static = true;
		mainElem = anotherMatrixDiag.mainElem;
		upElem = anotherMatrixDiag.upElem;
		downElem = anotherMatrixDiag.downElem;
		mainDiag = upDiag = downDiag = Vec();
	}
	else {
		mainElem = upElem = downElem = 0;
		mainDiag = anotherMatrixDiag.mainDiag;
		upDiag = anotherMatrixDiag.upDiag;
		downDiag = anotherMatrixDiag.downDiag;
	}

	if (anotherMatrixDiag.is_full == true) {
		is_full = true;
		matrix = (double **)malloc(sizeof(double *) * size);
		for (int i = 0; i < size; i++)
			matrix[i] = (double *)malloc(sizeof(double) * size);
		matrix[0][0] = anotherMatrixDiag.matrix[0][0];
		matrix[0][1] = anotherMatrixDiag.matrix[0][1];
		for (int i = 2; i < size; i++)
			matrix[0][i] = 0.0;
		int tmp = size - 1, tmp2 = size - 2;
		matrix[tmp][tmp] = anotherMatrixDiag.matrix[tmp][tmp];
		matrix[tmp][tmp2] = anotherMatrixDiag.matrix[tmp][tmp2];
		for (int i = 0; i < tmp2; i++)
			matrix[tmp][i] = 0.0;
		for (int i = 1; i < tmp; i++)
			for (int j = 0; j < size; j++)
				matrix[i][j] = 0.0;
		for (int i = 1; i < tmp; i++) {
			matrix[i][i] = anotherMatrixDiag.matrix[i][i];
			matrix[i][i - 1] = anotherMatrixDiag.matrix[i][i - 1];
			matrix[i][i + 1] = anotherMatrixDiag.matrix[i][i + 1];
		}
	}
	else {
		matrix = 0;
	}

}

MatrixDiag & MatrixDiag::operator=(const MatrixDiag & anotherMatrixDiag) // CHECK & CORRECT!
{
	if (&anotherMatrixDiag == this)
		return *this;

	if (anotherMatrixDiag.is_static) {
		is_static = true;
		mainDiag = upDiag = downDiag = Vec();
		mainElem = anotherMatrixDiag.mainElem;
		upElem = anotherMatrixDiag.upElem;
		downElem = anotherMatrixDiag.downElem;
	}
	else {
		is_static = false;
		mainElem = upElem = downElem = 0;
		mainDiag = anotherMatrixDiag.mainDiag;
		upDiag = anotherMatrixDiag.upDiag;
		downDiag = anotherMatrixDiag.downDiag;
	}
	if (anotherMatrixDiag.is_full) {
		is_full = true;

		if (size != anotherMatrixDiag.size) {
			if (size != 0) {
				for (int i = 0; i < size; i++)
					free(matrix[i]);
				free(matrix);
			}
			size = anotherMatrixDiag.size;
			matrix = (double **)malloc(sizeof(double *) * size);
			for (int i = 0; i < size; i++)
				matrix[i] = (double *)malloc(sizeof(double) * size);
		}
		
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				matrix[i][j] = anotherMatrixDiag.matrix[i][j];
	}
	else {
		if (is_full && (size != 0)) {
			for (int i = 0; i < size; i++)
				free(matrix[i]);
			free(matrix);
			
		}
		is_full = false;
		size = anotherMatrixDiag.size;
	}

	return *this;
}

MatrixDiag::~MatrixDiag()
{
	if (size == 0) {
		return;
	}

	if (is_full) {
		for (int i = 0; i < size; i++)
			free(matrix[i]);
		free(matrix);
	}
}

bool MatrixDiag::isEmpty()
{
	return (size == 0);
}

Vec MatrixDiag::sweep(Vec F)
{
	Vec X(size);
	Vec alfa(size);
	Vec beta(size);
	//__int64 ctr1 = 0, ctr2 = 0, freq = 0;
	//QueryPerformanceCounter((LARGE_INTEGER *)&ctr1);
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
	//QueryPerformanceCounter((LARGE_INTEGER *)&ctr2);
	//QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
	//double wtime((ctr2 - ctr1) * 1.0 / freq);
	//ofstream fout("wtime1.txt", std::fstream::app);
	//fout << wtime << endl;
	//fout.close();

	return X;
}
/*
double *MatrixDiag::sweepOMP(double *F, int nOfThreads)
{
	int curNOfThreads = omp_get_num_threads();
	omp_set_num_threads(nOfThreads);
	double *X = (double *)malloc(sizeof(double) * size);
	for (int i = 0; i < size; i++)
		X[i] = 0.0;

	int m = size / nOfThreads;
#pragma omp parallel
	{
		 int k = omp_get_thread_num();
		 int start = k * m + 1, end = (k + 1) * m - 1;
		 int tmp = 0;
		 double koef = 0.0;
		 
		 if (k == 0) { //  для верхнего квадрата
			for (int i = start; i <= end; i++) {
				koef = downDiag[i] / mainDiag[i - 1];
				matrix[i][i - 1] = 0;
				matrix[i][i] -= koef;
				F[i] -= F[i - 1] * koef;
			}
//#pragma omp barrier
			start = (k + 1) * m - 3;
			end = 0;
			tmp = start + 2;
			for (int i = start; i >= end; i--) {
				koef = upDiag[i] / matrix[i + 1][i + 1];
				matrix[i][i + 1] = 0;
				matrix[i][tmp] = -matrix[i + 1][tmp] * koef;
				F[i] -= F[i + 1] * koef;
			}
//#pragma omp barrier
			start = k * m;
			end = (k + 1) * m - 1;
			for (int i = start; i < end; i++) {
				X[i] = (F[i] - matrix[i][end]) / matrix[i][i];
			}
		 }
		 else if (k == nOfThreads - 1) { //  для нижнего квадрата
			 tmp = start - 2;
			 for (int i = start; i <= end; i++) {
				 koef = downDiag[i] / mainDiag[i - 1];
				 matrix[i][i - 1] = 0;
				 matrix[i][i] -= koef;
				 matrix[i][tmp] = -matrix[i - 1][tmp] * koef;
				 F[i] -= F[i - 1] * koef;
			 }
//#pragma omp barrier
			 start = (k + 1) * m - 3;
			 end = k * m - 1;
			 tmp = start + 2;
			 for (int i = start; i >= end; i--) {
				 koef = upDiag[i] / matrix[i + 1][i + 1];
				 matrix[i][i + 1] = 0;
				 matrix[i][tmp] = -matrix[i + 1][tmp] * koef;
				 matrix[i][end] -= -matrix[i + 1][end] * koef;
				 F[i] -= F[i + 1] * koef;
			 }
//#pragma omp barrier
			 start = k * m;
			 end = (k + 1) * m - 1;
			 for (int i = start; i < end; i++) {
				 X[i] = (F[i] - matrix[i][end] - matrix[i][start - 1]) / matrix[i][i];
			 }
		 }
		 else { //  для внутренних квадратов
			 tmp = start - 2;
			 for (int i = start; i <= end; i++) {
				 koef = downDiag[i] / mainDiag[i - 1];
				 matrix[i][i - 1] = 0;
				 matrix[i][i] -= koef;
				 matrix[i][tmp] = -matrix[i - 1][tmp] * koef;
				 F[i] -= F[i - 1] * koef;
			 }
//#pragma omp barrier
			 start = (k + 1) * m - 3;
			 end = k * m - 1;
			 tmp = start + 2;
			 for (int i = start; i >= end; i--) {
				 koef = upDiag[i] / matrix[i + 1][i + 1];
				 matrix[i][i + 1] = 0;
				 matrix[i][tmp] = -matrix[i + 1][tmp] * koef;
				 matrix[i][end] -= -matrix[i + 1][end] * koef;
				 F[i] -= F[i + 1] * koef;
			 }
//#pragma omp barrier
			 start = k * m;
			 end = (k + 1) * m - 1;
			 for (int i = start; i < end; i++) {
				 X[i] = (F[i] - matrix[i][end] - matrix[i][start - 1]) / matrix[i][i];
			 }
		 }
	 
	}

	vector<double> uptmp(nOfThreads), downtmp(nOfThreads), maintmp(nOfThreads), ftmp(nOfThreads);
	maintmp[0] = matrix[m - 1][m - 1];
	maintmp[nOfThreads - 1] = matrix[size - 1][size - 1];
	upDiag[0] = matrix[m - 1][2 * m - 1];
	downtmp[nOfThreads - 1] = matrix[size - 1][size - 1 - m];
	for (int i = 1; i < nOfThreads - 1; i++) {
		maintmp[i] = matrix[(i + 1) * m - 1][(i + 1) * m - 1];
		upDiag[i] = matrix[(i + 1) * m - 1][(i + 2) * m - 1];
		downtmp[i] = matrix[(i + 1) * m - 1][i * m - 1];
	}
	for (int i = 0; i < nOfThreads; i++) {
		ftmp[i] = F[(i + 1) * m - 1];
	}
	MatrixDiag L(nOfThreads, maintmp, uptmp, downtmp);
	vector<double> xtmp = L.sweep(ftmp);
	for (int i = 0; i < nOfThreads; i++) {
		X[(i + 1) * m - 1] = xtmp[i];
	}

	omp_set_num_threads(curNOfThreads);
	return X;
}

double *MatrixDiag::sweepOMP1(double *F)
{
	//int curNOfThreads = omp_get_num_threads();
	//omp_set_num_threads(2);
	//double wtime0 = 0.0, wtime1 = 0.0;
	vector<double> X(size, 0.0);
	vector<double> alfa(size, 0.0);
	vector<double> beta(size, 0.0);
	vector<double> ksi(size, 0.0);
	vector<double> eta(size, 0.0);
	int p = size / 2;

#pragma omp parallel shared(X, alfa, beta, ksi, eta, p) num_threads(2)
{
#pragma omp sections //firstprivate(wtime0, wtime1)
	{
#pragma omp section
		{
			//  thread 1
			//wtime0 = omp_get_wtime();
			int itmp = 0, ptmp = p + 1;
			alfa[1] = -upDiag[0] / mainDiag[0];
			beta[1] = F[0] / mainDiag[0];
			for (int i = 2; i <= ptmp; i++) {
				itmp = i - 1;
				alfa[i] = -upDiag[itmp] / (downDiag[itmp] * alfa[itmp] + mainDiag[itmp]);
				beta[i] = (F[itmp] - downDiag[itmp] * beta[itmp]) / (downDiag[itmp] * alfa[itmp] + mainDiag[itmp]);
			}
			//wtime0 = omp_get_wtime() - wtime0;
		}
#pragma omp section
		{
			//  thread 2
			//wtime0 = omp_get_wtime();
			int tmp = size - 1;
			int itmp = 0;
			ksi[tmp] = -downDiag[tmp] / mainDiag[tmp];
			eta[tmp] = F[tmp] / mainDiag[tmp];
			for (int i = tmp - 1; i > p; i--) {
				itmp = i + 1;
				ksi[i] = -downDiag[i] / (upDiag[i] * ksi[itmp] + mainDiag[i]);
				eta[i] = (F[i] - upDiag[i] * eta[itmp]) / (upDiag[i] * ksi[itmp] + mainDiag[i]);
			}
			//wtime0 = omp_get_wtime() - wtime0;
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
			//wtime1 = omp_get_wtime();
			int itmp = 0;
			for (int i = p - 1; i >= 0; i--) {
				X[i] = alfa[itmp] * X[itmp] + beta[itmp];
			}
			//wtime1 = omp_get_wtime() - wtime1;
		}
#pragma omp section
		{
			//  thread 2
			//wtime1 = omp_get_wtime();
			for (int i = p + 1; i < size; i++) {
				X[i] = ksi[i] * X[i - 1] + eta[i];
			}
			//wtime1 = omp_get_wtime() - wtime1;
		}
	}
}
	//omp_set_num_threads(curNOfThreads);
	//ofstream fout("wtime2.txt", std::fstream::app);
	//fout << wtime0 + wtime1 << endl;
	//fout.close();
	return X;
}

*/
