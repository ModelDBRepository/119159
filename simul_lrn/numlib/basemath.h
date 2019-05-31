#ifndef BASEMATH_H
#define BASEMATH_H

#include <cassert>
#include <cmath>
#define epsilon 1e-20


template <class Type>
inline Type max(Type, Type);

template <class Type>
inline Type min(Type, Type);

template <class Type>
inline Type abs(Type);


template <class Type>
inline Type ipower(Type, int);

double binomial(double, double);
long double binomial(long double, long double);
int binomial(int, int);

double p_binomial(double, double);
long double p_binomial(long double, long double);

#define sgn(X) (X>0? 1:-1)

#include "basemath.cc"
#endif
