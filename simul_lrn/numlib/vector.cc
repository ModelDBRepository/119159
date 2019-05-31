

#include "vector.h"

template <class Type>
vector<Type>::vector()
{
  n= 0;
  data= NULL;
}

template <class Type>
vector<Type>::vector(int sz)
{
  assert(sz >= 0);
  n= sz;
  if (sz == 0)
  {
    data= NULL;
  }
  else
  {
    data= new Type[n];
    for(int i= 0; i < n; i++)
    {
      data[i]= 0;
    }  
  }
}

template <class Type>
vector<Type>::vector(const vector<Type>& vec)
{
  n= vec.n;
  if (n == 0)
  {
    data= NULL;
  }
  else
  {
    data= new Type[n];
    for(int i= 0; i < n; i++)
    {
      data[i]= vec.data[i];
    }
  }
}

template <class Type>
vector<Type>::vector(int sz, Type *dat)
{
  assert(sz >= 0);
  n= sz;
  if (n == 0)
  {
    data= NULL;
  }
  else
  {
    data= new Type[n];
    for(int i= 0; i < n; i++)
    {
      data[i]= dat[i];
    }
  }
}

template <class Type>
vector<Type>::vector(int sz, char *dat)
{
  assert(sz >= 0);
  n= sz;
  if (n == 0) data= NULL;
  else
  {
    data= new Type[n];
    istrstream inp(dat);
    for(int i= 0; i < n; i++) inp >> data[i];
  }
}

template <class Type>
vector<Type>::~vector()
{
  if (n != 0)
  {
    delete[] data;
  }
}

template <class Type>
void vector<Type>::resize(int sz)
{
  assert(sz >= 0);
  if (n != 0)
  {
    delete[] data;
  }
  n= sz;
  if (n == 0)
  {
    data= NULL;
  }
  else
  {
    data= new Type[n];
    for(int i= 0; i < n; i++)
    {
      data[i]= 0;
    }
  }
}

template <class Type>
int vector<Type>::dim()
{
  return n;
}

template <class Type>
void vector<Type>::set(int i, Type dat)
{
  assert(i < n);
  data[i]= dat;
}

template <class Type>
const vector<Type>& vector<Type>::operator=(const vector<Type>& vec)
{
  if (n != vec.n)
  {
    if ( n != 0)
    {
      delete[] data;
    }
    n= vec.n;
    if (n == 0)
    {
      data= NULL;
    }
    else
    {
      data= new Type[n];
    }
  }
  for(int i= 0; i < n; i++)
  {
    data[i]= vec.data[i];
  }
  return *this;
}


template <class Type>
const vector<Type>& vector<Type>::operator=(const char *dat)
{
  if (n != 0) 
  {
    istrstream inp(dat);
    for(int i= 0; i < n; i++)
    {
      inp >> data[i];
    }
  }
  return *this;
}


template <class Type>
const vector<Type>& vector<Type>::operator+=(const vector<Type>& vec)
{
  for(int i= 0; i < n; i++)
  {
    data[i]+= vec.data[i];
  }
  return *this;
}


template <class Type>
const vector<Type>& vector<Type>::operator-=(const vector<Type>& vec)
{
  for(int i= 0; i < n; i++)
  {
    data[i]-= vec.data[i];
  }
  return *this;
}

template <class Type>
Type vector<Type>::operator[](int i)
{
  if (i < n)
  { 
    return data[i];
  }
  else
  {
    return 0;
  }
}

template <class Type>
int vector<Type>::operator==(const vector<Type>& b)
{
  if (n != b.n)
  {
    return 0;
  }
  for(int i= 0; i < n; i++)
  {
    if (data[i] != b.data[i])
    {
      return 0;
    }
  }
  return 1;
}

template <class Type>
int vector<Type>::approx_equal(const vector<Type>& b)
{
  if (n != b.n)
  {
    return 0;
  }
  for(int i= 0; i < n; i++)
  {
    if (abs(data[i]-b.data[i]) > epsilon)
    {
      return 0;
    }
  }
  return 1;
}


template <class Type>
int vector<Type>::operator!=(const vector<Type>& b)
{
  return (!(*this == b));
}


template <class Type>
vector<Type> operator+(const vector<Type>& a, const vector<Type>& b)
{
  assert(a.n == b.n);

  vector<Type> c(a);
  c+= b;
  return c;
}
  
template <class Type>
vector<Type> operator-(const vector<Type>& a, const vector<Type>& b)
{
  assert(a.n == b.n);
  
  vector<Type> c(a);
  c-= b;
  return c;
}

template <class Type>
Type operator*(const vector<Type>& a, const vector<Type>& b)
{
  assert(a.n == b.n);

  Type res= 0;
  for (int i= 0; i < a.n; i++)
  {
    res+= a.data[i]*b.data[i];
  }
  return res;
}


template <class Type>
vector<Type> operator*(Type a, const vector<Type>& b)
{
  vector<Type> c(b);
  for (int i= 0; i < b.n; i++)
  {
    c.data[i]*= a;
  }
  return c;
}


template <class Type>
istream& operator>>(istream& is, vector<Type>& a)
{
  Type dat;
  for(int i= 0; i < a.n; i++)
  {
    is >> dat;
    a.data[i]= dat;
  }
  return is;
}

template <class Type>
ostream& operator<<(ostream& os, const vector<Type>& a)
{
  for(int i= 0; i < a.n; i++)
  {
    os << a.data[i] << " ";
  }
  return os;
}


template <class Type>
Type norm1(vector<Type> a)
{
  assert(a.n > 0);
  
  Type sum= abs(a.data[0]);
  for (int i= 1; i < a.n; i++)
  {
    sum+= abs(a.data[i]);
  }
  return sum;
}


template <class Type>
Type norm2(vector<Type> a)
{
  assert(a.n > 0);
  
  Type sum= a.data[0]*a.data[0];
  for (int i= 1; i < a.n; i++)
  {
    sum+= a.data[i]*a.data[i];
  }
  sum= (Type) sqrt(sum);
  return sum;
}


