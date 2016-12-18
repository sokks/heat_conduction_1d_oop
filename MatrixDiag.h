#pragma once
#include <vector>
#include <cmath>
#include <omp.h>
#include <Windows.h>
#include <fstream>
#include <iostream>

using namespace std;

class MatrixDiag
{
	int size;
	double *mainDiag, *upDiag, *downDiag;
	bool is_static = false;
	double mainElem, upElem, downElem;
	bool is_full = false;
	double **matrix;
public:
	MatrixDiag();
	MatrixDiag(int n, double diagElement, double upElement, double downElement, bool stat = false, bool full = false);
	MatrixDiag(int n, double *diag, double *up, double *down);
	MatrixDiag(const MatrixDiag & anotherMatrixDiag);
	MatrixDiag & operator= (const MatrixDiag & anotherMatrixDiag);
	~MatrixDiag();

	bool isEmpty();
	double *sweep(double *F);
	double *sweepOMP(double *F, int nOfThreads = 4);
	double *sweepOMP1(double *F);
};

