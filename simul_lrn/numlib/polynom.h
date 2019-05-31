
#ifndef POLYNOM_H
#define POLYNOM_H

#include <iostream>
#include <assert.h>
#include <String.h>
#include <math.h>

#include "TNmath.h"

template <class Type>
class polynomial
{
    friend istream& operator>>(istream&, polynomial<Type>&);                  
    friend ostream& operator<<(ostream&, polynomial<Type>&);                  

private:
    int dimen, ord;
    tnvector<Type> coeff;
    tnvector<Type> org;

public:
    polynomial();
    polynomial(int, int);
    polynomial(int, int, tnvector<Type>&, tnvector<Type>&);
    polynomial(polynomial<Type>&);
    polynomial(istream&);    
    ~polynomial();
    void resize(int, int);
    int dim();
    int order();
    int coeff_size();
    Type value(tnvector<Type>);
    Type coefficient(tnvector<int>);
    void set_coeff_vec(tnvector<Type>);
    void set_coeff_direct(int pos, Type value);
    tnvector<Type> coeff_vec();
    void set_origin(tnvector<Type>);
    tnvector<Type> origin();
    void read(istream&);
    void nice_output(ostream&);
};

#include "polynom.cc"

#endif
