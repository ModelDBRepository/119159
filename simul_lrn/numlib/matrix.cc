

#include "matrix.h"

template <class Type>
matrix<Type>::matrix()
{
  m= 1;
  n= 1;
  data= new Type[1];
  data[0]= 0;
}

template <class Type>
matrix<Type>::matrix(int zsz, int ssz)
{
  m= zsz;
  n= ssz;
  data= new Type[m*n];
  for(int i= 0; i < m; i++)
  {
    for(int j= 0; j < n; j++)
    {
      data[i*n+j]= 0;
    }
  }
}

template <class Type>
matrix<Type>::matrix(matrix<Type>& a)
{
  m= a.zdim();
  n= a.sdim();
  data= new Type[m*n];
  for(int i= 0; i < m; i++)
  {
    for(int j= 0; j < n; j++)
    {
      data[i*n+j]= a[i][j];
    }
  }
}


template <class Type>
matrix<Type>::matrix(int zsz, int ssz, Type *a)
{
  m= zsz;
  n= ssz;
  data= new Type[m*n];
  for(int i= 0; i < m; i++)
  {
    for(int j= 0; j < n; j++)
    {
      data[i*n+j]= a[i*n+j];
    }
  }
}

template <class Type>
matrix<Type>::~matrix()
{
  delete[] data;
}


template <class Type>
void matrix<Type>::resize(int zsz, int ssz)
{
  m= zsz;
  n= ssz;
  delete[] data;
  data= new Type[m*n];
  for(int i= 0; i < m; i++)
  {
    for(int j= 0; j < n; j++)
    {
      data[i*n+j]= 0;
    }
  }
}  


template <class Type>
int matrix<Type>::zdim()
{
  return m;
}
template <class Type>

int matrix<Type>::sdim()
{
  return n;
}

template <class Type>
void matrix<Type>::set(int i, int j, Type dat)
{
  assert(i < m);
  assert(j < n);

  data[i*n+j]= dat;
}

template <class Type>
void matrix<Type>::swap_lines(int i, int k)
{
  assert(i < m);
  assert(k < m);

  Type buf= 0;
  for(int j= 0; j < n; j++)
  {
    buf= data[i*n+j];
    data[i*n+j]= data[k*n+j];
    data[k*n+j]= buf;
  }
}


template <class Type>
void matrix<Type>::swap_cols(int j, int k)
{
  assert(j < n);
  assert(k < n);

  Type buf= 0;
  for(int i= 0; i < m; i++)
  {
    buf= data[i*n+j];
    data[i*n+j]= data[i*n+k];
    data[i*n+k]= buf;
  }
}


template <class Type>
int matrix<Type>::invert(matrix<Type>& b)
{
  Type fac;
  assert(m == n);

  matrix<Type> a(*this);
  
  for(int i= 0; i < m; i++)
  {
    for(int j= 0; j < n; j++)
    {
      if (i == j)
      {
	b.set(i,j,1);
      }
      else
      {
	b.set(i,j,0);
      }
    }
  }
  
  for(int i= 0; i < m; i++)
  {
    if (fabs(a[i][i]) < epsilon)
    {
      int l= i;
      while((l < m) && (fabs(a[l][i]) < epsilon))
      {
	l++;
      }
      if (l == m)
      {
	return 0;
      }
      else
      {
	a.swap_lines(i,l);
	b.swap_lines(i,l);
      }
    }
    for(int k= i; k < m; k++)
    {
      if (fabs(a[k][i]) > epsilon)
      {
	fac= 1/a[k][i];
	for(int j= 0; j < n; j++)
	{
	  a.set(k,j,a[k][j]*fac);
	  b.set(k,j,b[k][j]*fac);
	}
      }
    }
    for(int k= i+1; k < m; k++)
    {
      if (fabs(a[k][i]) > epsilon)
      {
	for(int j= 0; j < n; j++)
	{
	  a.set(k,j,a[k][j]-a[i][j]);
	  b.set(k,j,b[k][j]-b[i][j]);
	}
      }
    }
  }
  for(int i= m-1; i > 0; i--)
  {
    for(int k= 0; k < i; k++)
    {
      if (fabs(a[k][i]) > epsilon)
      {
	fac= a[k][i];
	for(int j= 0; j < n; j++)
	{
	  a.set(k,j,a[k][j]-a[i][j]*fac);
	  b.set(k,j,b[k][j]-b[i][j]*fac);
	}
      }
    }
  }
  return 1;
}



template <class Type>
const matrix<Type>& matrix<Type>::operator=(const matrix<Type>& a)
{
  if ((m != a.m) || (n != a.n))
  {
    m= a.m;
    n= a.n;
    delete[] data;
    data= new Type[m*n];
  }
  for(int i= 0; i < m; i++)
  {
    for(int j= 0; j < n; j++)
    {
      data[i*n+j]= a.data[i*n+j];
    }
  }
  return *this;
}

template <class Type>
const matrix<Type>& matrix<Type>::operator+=(const matrix<Type>& a)
{
  for(int i= 0; i < m; i++)
  {  
    for(int j= 0; j < n; j++)
    {
      data[i*n+j]+= a.data[i*n+j];
    }
  }
  return *this;
}


template <class Type>
const matrix<Type>& matrix<Type>::operator-=(const matrix<Type>& a)
{
  for(int i= 0; i < m; i++)
  {
    for(int j= 0; j < n; j++)
    {
      data[i*n+j]-= a.data[i*n+j];
    }
  }
  return *this;
}

template <class Type>
tnvector<Type> matrix<Type>::operator[](int i)
{
  if (i < m)
  { 
    tnvector<Type> vec(n,&data[i*n]);
    return vec;
  }
  else
  {
    tnvector<Type> vec(n);    
    return vec;
  }
}


template <class Type>
int matrix<Type>::operator==(const matrix<Type>& b)
{
  if ((m != b.m) || (n != b.n))
  {
    return 0;
  }
  for(int i= 0; i < m; i++)
  {
    for(int j= 0; j < n; j++)
    {
      if (data[i*n+j] != b.data[i*n+j])
      {
	return 0;
      }
    }
  }
  return 1;
}

template <class Type>
int matrix<Type>::operator!=(const matrix<Type>& b)
{
  return (!(*this == b));
}


template <class Type>
matrix<Type> operator+(const matrix<Type>& a, const matrix<Type>& b)
{
  assert((a.m == b.m) && (a.n == b.n));

  matrix<Type> c(a);
  c+= b;
  return c;
}
  
template <class Type>
matrix<Type> operator-(const matrix<Type>& a, const matrix<Type>& b)
{
  assert((a.m == b.m) && (a.n == b.n));

  matrix<Type> c(a);
  c-= b;
  return c;
}

template <class Type>
matrix<Type> operator*(const matrix<Type>& a, const matrix<Type>& b)
{
  assert(a.n == b.m);

  matrix<Type> c(a.m, b.n);
  Type buf;
  for (int i= 0; i < a.m; i++)
  {
    for (int j= 0; j < b.n; j++)
    {
      buf= 0;
      for (int k= 0; k < a.n; k++)
      {
	buf+= a.data[i*a.n+k]*b.data[k*b.n+j];
      }
      c.data[i*c.n+j]= buf;
    }
  }
  return c;
}


template <class Type>
tnvector<Type> operator*(const matrix<Type>& a, const tnvector<Type>& b)
{
  assert(a.n == b.n);

  tnvector<Type> c(a.m);
  Type buf;

  for(int i= 0; i < a.m; i++)
  {
    buf= 0;
    for(int j= 0; j < a.n; j++)
    {
      buf+= a.data[i*a.n+j]*b.data[j];
    } 
    c.data[i]= buf;
  }
  return c;
}


template <class Type>
matrix<Type> operator*(Type a, const matrix<Type>& b)
{
  matrix<Type> c(b);
  for (int i= 0; i < b.m; i++)
  {
    for (int j= 0; j < b.n; j++)
    {
      c.data[i*c.n+j]*= a;
    }
  }
  return c;
}

template <class Type>
istream& operator>>(istream& is, matrix<Type>& a)
{
  Type dat;
  for(int i= 0; i < a.zdim(); i++)
  {
    for(int j= 0; j < a.sdim(); j++)
    {
      is >> dat;
      a.set(i, j, dat);
    }
  }
  return is;
}

template <class Type>
ostream& operator<<(ostream& os, matrix<Type> a)
{
  for(int i= 0; i < a.zdim(); i++)
  {
    for(int j= 0; j < a.sdim(); j++)
    {
      os << a[i][j] << " ";
    }
    os << endl;
  }
  return os;
}

