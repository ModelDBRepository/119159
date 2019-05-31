//--------------------------------------------------------------------------
// Author: Thomas Nowotny
//
// Institute: Institut fuer Theoretische Physik
//            Augustusplatz 10-11
//            04109 Leipzig
//
// email to:  nowotny@itp.uni-leipzig.de
//
// initial version: 8/99
// last change: 09/00
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
// the type eld was defined to implement simple interval arithmetic for
// long double variables.
// note that the definition of <, >, <=, >=, ==, != are somewhat arbitrary
// and that the choices have been made in quite different ways for each.
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
// some usage information:
// In order to use the extended long doubles (eld) for the base number type
// type the constants null (the number 0), one (the number 1), __eld_eps
// (the smallest number with one+__eld_eps != one) and __eld_pi (the number
// pi in the representation of type) must be defined globally
//--------------------------------------------------------------------------

#ifndef ELD_H
#define ELD_H

#include <iostream>
#include <math.h>
#include "mycmath"
#include "basemath.h"

template <class type>
class eld
{
 public:
    type xm, x, xp;
    eld();
    eld(type);
    eld(type,type,type);
    eld(const eld&);
    inline eld& operator=(const eld);
    inline eld& operator=(const type);
    ~eld() { }
    inline eld inv();
    inline eld& operator+=(const eld);
    inline eld& operator-=(const eld);
    inline eld& operator*=(const eld);
    inline eld& operator/=(const eld);
    inline eld& operator+=(const type);
    inline eld& operator-=(const type);
    inline eld& operator*=(const type);
    inline eld& operator/=(const type);
};

template <class type>
inline eld<type> operator-(const eld<type>);
template <class type>
inline eld<type> operator+(const eld<type>, const eld<type>);
template <class type>
inline eld<type> operator-(const eld<type>, const eld<type>);
template <class type>
inline eld<type> operator*(const eld<type>, const eld<type>);
template <class type>
inline eld<type> operator/(const eld<type>, const eld<type>);

template <class type>
inline eld<type> operator+(const type, const eld<type>);
template <class type>
inline eld<type> operator-(const type, const eld<type>);
template <class type>
inline eld<type> operator*(const type, const eld<type>);
template <class type>
inline eld<type> operator/(const type, const eld<type>);

template <class type>
inline eld<type> operator+(const eld<type>, const type);
template <class type>
inline eld<type> operator-(const eld<type>, const type);
template <class type>
inline eld<type> operator*(const eld<type>, const type);
template <class type>
inline eld<type> operator/(const eld<type>, const type);

template <class type>
inline eld<type> sqrt(const eld<type>);
template <class type>
inline eld<type> sinh(const eld<type>);
template <class type>
inline eld<type> cosh(const eld<type>);
template <class type>
inline eld<type> tanh(const eld<type>);
template <class type>
inline eld<type> log(const eld<type>);
template <class type>
inline eld<type> exp(const eld<type>);
template <class type>
inline eld<type> atan(const eld<type>);
/*
template <class type>
inline eld<type> cos(const eld<type>);
template <class type>
inline eld<type> sin(const eld<type>);
*/
template <class type>
inline eld<type> tan(const eld<type>);
template <class type>
inline eld<type> abs(const eld<type>);
template <class type>
inline eld<type> pow(const eld<type>, const eld<type>);
template <class type>
inline eld<type> pow(const eld<type>, const type);
template <class type>
inline eld<type> pow(const eld<type>, int);
template <class type>
inline eld<type> pow(const type, const eld<type>);
template <class type>
inline eld<type> atanh(const eld<type>);
template <class type>
inline eld<type> floor(const eld<type>);

template <class type>
inline int operator==(const eld<type>, const eld<type>);
template <class type>
inline int operator!=(const eld<type>, const eld<type>);
template <class type>
inline int operator<(const eld<type>, const eld<type>);
template <class type>
inline int operator<=(const eld<type>, const eld<type>);
template <class type>
inline int operator>(const eld<type>, const eld<type>);
template <class type>
inline int operator>=(const eld<type>, const eld<type>);

template <class type>
inline int operator<<(const eld<type>, const eld<type>);

template <class type>
inline int operator==(const eld<type>, const type);
template <class type>
inline int operator!=(const eld<type>, const type);
template <class type>
inline int operator<(const eld<type>, const type);
template <class type>
inline int operator<=(const eld<type>, const type);
template <class type>
inline int operator>(const eld<type>, const type);
template <class type>
inline int operator>=(const eld<type>, const type);

template <class type>
inline int operator==(const type, const eld<type>);
template <class type>
inline int operator!=(const type, const eld<type>);
template <class type>
inline int operator<(const type, const eld<type>);
template <class type>
inline int operator<=(const type, const eld<type>);
template <class type>
inline int operator>(const type, const eld<type>);
template <class type>
inline int operator>=(const type, const eld<type>);

template <class type>
inline eld<type> min(const eld<type>, const type);
template <class type>
inline eld<type> min(const type, const eld<type>);
template <class type>
inline eld<type> max(const eld<type>, const type);
template <class type>
inline eld<type> max(const type, const eld<type>);

template <class type>
inline long conv_to_long(const eld<type>);
template <class type>
inline short conv_to_short(const eld<type>);

template <class type>
inline eld<type> trunc(const eld<type>);

#include "eld.cc"

#endif





