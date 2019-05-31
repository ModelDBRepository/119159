#ifndef TNVECTOR_H
#define TNVECTOR_H

#include <iostream>
#include <sstream>
#include <assert.h>

template <class Type>
class tnvector;                             // forward declaration

#include "matrix.h"

template <class Type>
class tnvector;

template <class Type>
tnvector<Type> operator+(const tnvector<Type>&, const tnvector<Type>&);

template <class Type>
tnvector<Type> operator-(const tnvector<Type>&, const tnvector<Type>&);

template <class Type>
tnvector<Type> operator*(const matrix<Type>&, const tnvector<Type>&);

template <class Type>
Type operator*(const tnvector<Type>&, const tnvector<Type>&);

template <class Type>
tnvector<Type> operator*(Type, const tnvector<Type>&);

template <class Type>
Type norm1(tnvector<Type>);

template <class Type>
Type norm2(tnvector<Type>);

template <class Type>
istream& operator>>(istream&, tnvector<Type>&);

template <class Type>
ostream& operator<<(ostream&, const tnvector<Type>&);    

template <class Type>
class tnvector
{
  friend tnvector<Type> operator+<>(const tnvector<Type>&, const tnvector<Type>&);
  friend tnvector<Type> operator-<>(const tnvector<Type>&, const tnvector<Type>&);
  friend tnvector<Type> operator*<>(const matrix<Type>&, const tnvector<Type>&);
  friend Type operator*<>(const tnvector<Type>&, const tnvector<Type>&);
  friend tnvector<Type> operator*<>(Type, const tnvector<Type>&);

  friend Type norm1<>(tnvector<Type>);
  friend Type norm2<>(tnvector<Type>);

  friend istream& operator>><>(istream&, tnvector<Type>&);
  friend ostream& operator<<<>(ostream&, const tnvector<Type>&);    
  
protected:
    int n;
    Type *data;
    
public:
    tnvector();
    tnvector(int);
    tnvector(const tnvector<Type>&);
    tnvector(int, Type*);
    tnvector(int, char*);
    ~tnvector();
    void resize(int);
    int dim();
    void set(int, Type);
    const tnvector<Type>& operator=(const tnvector<Type>&);
    const tnvector<Type>& operator=(const char *dat);
    const tnvector<Type>& operator+=(const tnvector<Type>&);
    const tnvector<Type>& operator-=(const tnvector<Type>&);  
    Type operator[](int i);
    int operator==(const tnvector<Type>&);
    int operator!=(const tnvector<Type>&);
    int approx_equal(const tnvector<Type>&);
};

template <class Type>
tnvector<Type> operator+(const tnvector<Type>&, const tnvector<Type>&);

template <class Type>
tnvector<Type> operator-(const tnvector<Type>&, const tnvector<Type>&);

template <class Type>
tnvector<Type> operator*(const matrix<Type>&, const tnvector<Type>&);

template <class Type>
Type operator*(const tnvector<Type>&, const tnvector<Type>&);

template <class Type>
tnvector<Type> operator*(Type, const tnvector<Type>&);

template <class Type>
Type norm1(tnvector<Type>);

template <class Type>
Type norm2(tnvector<Type>);

template <class Type>
istream& operator>>(istream&, tnvector<Type>&);

template <class Type>
ostream& operator<<(ostream&, const tnvector<Type>&);


#include "tnvector.cc"

#endif
