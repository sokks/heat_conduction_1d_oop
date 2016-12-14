#pragma once
#include "Vec.h"
#include <cmath>
#include <omp.h>
#include <Windows.h>
#include <fstream>
#include <iostream>
#include <cstdlib>

using namespace std;

class MatrixDiag
{
	int size;
	Vec mainDiag, upDiag, downDiag;
	bool is_static = false;
	double mainElem, upElem, downElem;
	bool is_full = false;
	double **matrix;
public:
	MatrixDiag();
	MatrixDiag(int n, double diagElement, double upElement, double downElement, bool stat = false, bool full = false);
	MatrixDiag(int n, Vec diag, Vec up, Vec down);
	MatrixDiag(const MatrixDiag & anotherMatrixDiag);
	MatrixDiag & operator= (const MatrixDiag & anotherMatrixDiag);
	~MatrixDiag();

	bool isEmpty();
	Vec sweep(Vec F);
	Vec sweepOMP(Vec F, int nOfThreads = 4) { return 0; }
	Vec sweepOMP1(Vec F) { return 0; }
};

