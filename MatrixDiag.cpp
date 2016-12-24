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

	return X;
}

vector<double> MatrixDiag::sweep(vector<double> F, int mod, double right, double alfa1, double beta1, 
															double koef, double k1, double k2)
{
	vector<double> X(size, 0.0);
	vector<double> alfa(size, 0.0);
	vector<double> beta(size, 0.0);

	alfa[1] = alfa1;
	beta[1] = beta1;

	for (int i = 2; i < size; i++) {
		alfa[i] = -upDiag[i - 1] / (downDiag[i - 1] * alfa[i - 1] + mainDiag[i - 1]);
		beta[i] = (F[i - 1] - downDiag[i - 1] * beta[i - 1]) / (downDiag[i - 1] * alfa[i - 1] + mainDiag[i - 1]);
	}
	
	if (mod == 1) {
		X[size - 1] = right;
	}
	else if (mod == 2) {
		X[size - 1] = (k1 + koef * beta[size - 1]) / (k2 - koef * alfa[size - 1]);
	}
	for (int i = size - 2; i >= 0; i--) {
		X[i] = alfa[i + 1] * X[i + 1] + beta[i + 1];
	}

	return X;
}
