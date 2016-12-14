#pragma once
#include <iostream>

using namespace std;

class Vec
{
	int len;
	double *v;
public:
	Vec(int l = 0, double x = 0.0);
	Vec(const Vec & anotherVec);
	Vec & operator= (const Vec & anotherVec);
	~Vec();

	double & operator[] (int i);

	friend ostream& operator<< (ostream & out, const Vec & myvec);
};

