
#ifndef DPOLYNOM_H
#define DPOLYNOM_H

#include <iostream>
#include <assert.h>
#include <String.h>
#include <math.h>

#include "TNmath.h"
#include "dclock.h"

template <class Type>
class dpolynom
{
    friend istream& operator>>(istream&, dpolynom<Type>&);                  
    friend ostream& operator<<(ostream&, dpolynom<Type>&);                  

private:
    int dimen, ord;
    tnvector<Type> coeff;
    tnvector<Type> org;

public:
    dpolynom();
    dpolynom(int, int);
    dpolynom(int, int, tnvector<Type>&, tnvector<Type>&);
    dpolynom(dpolynom<Type>&);
    dpolynom(istream&);    
    ~dpolynom();
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
    void nice_output(ostream& os);
};

#include "dpolynom.cc"

#endif
