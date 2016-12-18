#pragma once
#include <vector>
#include <cmath>
//#include <omp.h>
#include <Windows.h>
#include <fstream>

using namespace std;

class MatrixDiag
{
	int size;
	vector<double> mainDiag, upDiag, downDiag;
public:
	MatrixDiag();
	MatrixDiag(int n, double diagElement, double upElement, double downElement);
	MatrixDiag(const MatrixDiag & anotherMatrixDiag);
	MatrixDiag & operator= (const MatrixDiag & anotherMatrixDiag);
	~MatrixDiag();

	bool isEmpty();
	vector<double> sweep(vector<double> F);
	vector<double> sweep(vector<double> F, double left, double right);
};

