#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <assert.h>
#include <math.h>
#include "basemath3.h"

template <class Type>
class matrix;

template <class Type>
matrix<Type> operator+(const matrix<Type>&, const matrix<Type>&);

template <class Type>
matrix<Type> operator-(const matrix<Type>&, const matrix<Type>&);

template <class Type>
matrix<Type> operator*(const matrix<Type>&, const matrix<Type>&);

template <class Type>
tnvector<Type> operator*(const matrix<Type>&, const tnvector<Type>&);

template <class Type>
matrix<Type> operator*(Type, const matrix<Type>&);

template <class Type>
class matrix
{
    friend matrix<Type> operator+<>(const matrix<Type>&, const matrix<Type>&);
    friend matrix<Type> operator-<>(const matrix<Type>&, const matrix<Type>&);
    friend matrix<Type> operator*<>(const matrix<Type>&, const matrix<Type>&);
    friend tnvector<Type> operator*<>(const matrix<Type>&, const tnvector<Type>&);
    friend matrix<Type> operator*<>(Type, const matrix<Type>&);
    
protected:
    int m, n;
    Type *data;
    
public:
    matrix();
    matrix(int, int);
    matrix(matrix<Type>&);
    matrix(int, int, Type *);
    ~matrix();
    void resize(int, int);
    int zdim();
    int sdim();
    void set(int, int, Type);
    void swap_lines(int, int);
    void swap_cols(int,int);
    int invert(matrix<Type>&);
    const matrix<Type>& operator=(const matrix<Type>&);
    const matrix<Type>& operator+=(const matrix<Type>&);
    const matrix<Type>& operator-=(const matrix<Type>&);  
    tnvector<Type> operator[](int);
    int operator==(const matrix<Type>&);
    int operator!=(const matrix<Type>&);
};

#include "matrix.cc"

#endif
