#include "Vec.h"

Vec::Vec(int l, double x)
{
	len = l;
	v = (double *)malloc(sizeof(double) * len);
	for (int i = 0; i < len; i++)
		v[i] = x;
}

Vec::Vec(const Vec & anotherVec)
{
	len = anotherVec.len;
	v = (double *)malloc(sizeof(double) * len);
}

Vec & Vec::operator=(const Vec & anotherVec)
{
	if (this == &anotherVec)
		return *this;

	if (len != anotherVec.len) {
		if (len != 0)
			free(v);
		len = anotherVec.len;
		v = (double *)malloc(sizeof(double) * len);
	}

	for (int i = 0; i < len; i++)
		v[i] = anotherVec.v[i];

	return *this;
}

Vec::~Vec()
{
	if (len != 0)
		free(v);
}

double & Vec::operator[](int i)
{
	if (i >= len) {
		throw - 1;
	}
	return v[i];
}

ostream & operator<<(ostream & out, const Vec & myvec)
{
	for (int i = 0; i < myvec.len; i++)
		out << myvec.v[i] << " ";
	return out;
}
