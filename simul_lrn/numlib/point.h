#ifndef POINT_H
#define POINT_H

#include "TNmath.h"

template <class Type>
class point
{    
public:
    tnvector<Type> x;
    Type y;

    point();
    point(int);
    point(point<Type>&);
    ~point() { }

    int operator==(const point<Type>&);
};

template <class Type>
class simple_point
{    
public:
    Type x;
    Type y;

    simple_point();
    simple_point(simple_point<Type>&);
    ~simple_point() { }

    int operator==(const simple_point<Type>&);
};

#include "point.cc"

#endif
