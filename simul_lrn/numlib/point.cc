

#include "point.h"

template <class Type>
point<Type>::point()
{
  y= 0;
}

template <class Type>
point<Type>::point(int sz)
{
  x.resize(sz);
  y= 0;
}

template <class Type>
point<Type>::point(point<Type>& p)
{
  x= p.x;
  y= p.y;
}


template <class Type>
int point<Type>::operator==(const point<Type>& b)
{
  if ((x== b.x) && (y== b.y))
  {
    return 1;
  }
  else
  {
    return 0;
  }
}



template <class Type>
simple_point<Type>::simple_point()
{
  x= 0;
  y= 0;
}

template <class Type>
simple_point<Type>::simple_point(simple_point<Type>& p)
{
  x= p.x;
  y= p.y;
}


template <class Type>
int simple_point<Type>::operator==(const simple_point<Type>& b)
{
  if ((x== b.x) && (y== b.y))
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

template <class Type>
ostream& operator<<(ostream& os, const point<Type>& a)
{
  os << a.x;
  os << " " << a.y;

  return os;
}
