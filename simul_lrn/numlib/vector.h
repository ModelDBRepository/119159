#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <strstream>
#include <assert.h>

template <class Type>
class vector;                             // forward declaration

#include "matrix.h"

template <class Type>
class vector
{
  friend vector<Type> operator+<>(const vector<Type>&, const vector<Type>&);
  friend vector<Type> operator-<>(const vector<Type>&, const vector<Type>&);
  friend vector<Type> operator*<>(const matrix<Type>&, const vector<Type>&);
  friend Type operator*<>(const vector<Type>&, const vector<Type>&);
  friend vector<Type> operator*<>(Type, const vector<Type>&);
  friend istream& operator>><>(istream&, vector<Type>&);
  friend ostream& operator<<<>(ostream&, const vector<Type>&);    
  friend Type norm1<>(vector<Type>);
  friend Type norm2<>(vector<Type>);
  
protected:
    int n;
    Type *data;
    
public:
    vector();
    vector(int);
    vector(const vector<Type>&);
    vector(int, Type*);
    vector(int, char*);
    ~vector();
    void resize(int);
    int dim();
    void set(int, Type);
    const vector<Type>& operator=(const vector<Type>&);
    const vector<Type>& operator=(const char *dat);
    const vector<Type>& operator+=(const vector<Type>&);
    const vector<Type>& operator-=(const vector<Type>&);  
    Type operator[](int i);
    int operator==(const vector<Type>&);
    int operator!=(const vector<Type>&);
    int approx_equal(const vector<Type>&);
};

template <class Type>
vector<Type> operator+(const vector<Type>&, const vector<Type>&);

template <class Type>
vector<Type> operator-(const vector<Type>&, const vector<Type>&);

template <class Type>
vector<Type> operator*(const matrix<Type>&, const vector<Type>&);

template <class Type>
Type operator*(const vector<Type>&, const vector<Type>&);

template <class Type>
vector<Type> operator*(Type, const vector<Type>&);

template <class Type>
istream& operator>>(istream&, vector<Type>&);

template <class Type>
ostream& operator<<(ostream&, const vector<Type>&);

template <class Type>
Type norm1(vector<Type>);

template <class Type>
Type norm2(vector<Type>);

#include "vector.cc"

#endif
