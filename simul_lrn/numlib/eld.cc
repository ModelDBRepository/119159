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
// last change: 8/99
//--------------------------------------------------------------------------

#include "eld.h"

template <class type>
eld<type>::eld()
{
}

template <class type>
eld<type>::eld(type xi)
{
  xm= xi;
  x= xi;
  xp= xi;
  assert(xm <= xp);
}

template <class type>
eld<type>::eld(type xmi, type xi, type xpi)
{
  xm= xmi;
  x= xi;
  xp= xpi;
  assert(xm <= xp);
}

template <class type>
eld<type>::eld(const eld<type> &ei)
{
  xm= ei.xm;
  x= ei.x;
  xp= ei.xp;
  assert(xm <= xp);
}

template <class type>
inline eld<type>& eld<type>::operator=(const eld<type> ei)
{
  xm= ei.xm;
  x= ei.x;
  xp= ei.xp;
  assert(xm <= xp);
  return *this;
}

template <class type>
inline eld<type>& eld<type>::operator=(const type li)
{
  xm= li;
  x= li;
  xp= li;
  assert(xm <= xp);
  return *this;
}

template <class type>
inline eld<type>& eld<type>::operator+=(const eld<type> e2)
{
    xm+= e2.xm;
    xm-= abs(xm*__eld_eps);
    x+= e2.x;
    xp+= e2.xp;
    xp+= abs(xp*__eld_eps);
    assert(xm <= xp);
    return *this;
}

template <class type>
inline eld<type>& eld<type>::operator-=(const eld<type> e2)
{
  xm-= e2.xp;
  x-= e2.x;
  xp-= e2.xm;
  xm-= abs(xm*__eld_eps);
  xp+= abs(xp*__eld_eps);
  assert(xm <= xp);
  return *this;
}

template <class type>
inline eld<type>& eld<type>::operator*=(const eld<type> e2)
{
  type h1= xm*e2.xm;
  type h2= xm*e2.xp;
  type h3= xp*e2.xm;
  type h4= xp*e2.xp;
  
  xm= min(min(h1,h2),min(h3,h4));
  xm-= abs(xm*__eld_eps);
  x*= e2.x;
  xp= max(max(h1,h2),max(h3,h4));
  xp+= abs(xp*__eld_eps);
  assert(xm <= xp);
  return *this;
}

template <class type>
inline eld<type> eld<type>::inv()
{
  eld<type> ei;
  if ((xm <= 0) && (xp >= 0))
  {
    ei.xm= -__eld_max;
    ei.x= one;
    ei.xp= __eld_max;
  }
  else
  {
    ei.xm= one/xp;
    ei.xm-= abs(ei.xm*__eld_eps);
    ei.x= one/x;
    ei.xp= one/xm;
    ei.xp+= abs(ei.xp*__eld_eps);
  }
  assert(ei.xm <= ei.xp);
  return ei;
}
	
template <class type>
inline eld<type>& eld<type>::operator/=(const eld<type> e2)
{
  eld<type> en(e2);
  *this*= en.inv();
  assert(xm <= xp);
  return *this;
}

template <class type>
inline eld<type>& eld<type>::operator+=(const type l)
{
  xm+= l;
  xm-= abs(xm*__eld_eps);
  x+= l;
  xp+= l;
  xp+= abs(xp*__eld_eps);
  assert(xm <= xp);
  return *this;
}

template <class type>
inline eld<type>& eld<type>::operator-=(const type l)
{
  xm-= l;
  xm-= abs(xm*__eld_eps);
  x-= l;
  xp-= l;
  xp+= abs(xp*__eld_eps);
  assert(xm <= xp);
  return *this;
}

template <class type>
inline eld<type>& eld<type>::operator*=(const type l)
{
  if (l < 0)
  {
    type tmp= xm;
    xm= xp*l;
    xp= tmp*l;
  }
  else
  {
    xm*= l;
    xp*= l;
  }
  x*= l;
  xm-= abs(xm*__eld_eps);
  xp+= abs(xp*__eld_eps);
  assert(xm <= xp);
  return *this;
}

template <class type>
inline eld<type>& eld<type>::operator/=(const type l)
{
  if (l < 0)
  {
    type tmp= xm;
    xm= xp/l;
    xp= tmp/l;
  }
  else
  {
    xm/= l;
    xp/= l;
  }
  x/= l;
  xm-= abs(xm*__eld_eps);
  xp+= abs(xp*__eld_eps);
  assert(xm <= xp);
  return *this;
}


template <class type>
inline eld<type> operator-(const eld<type> e)
{
  eld<type> en;
  en.xm= -e.xp;
  en.x= -e.x;
  en.xp= -e.xm;
  
  assert(en.xm <= en.xp);
  return en;
}

template <class type>
inline eld<type> operator+(const eld<type> e1, const eld<type> e2)
{
  eld<type> en(e1);
  en+= e2;

  assert(en.xm <= en.xp);
  return en;
}

template <class type>
inline eld<type> operator-(const eld<type> e1, const eld<type> e2)
{
  eld<type> en(e1);
  en-= e2;
  
  assert(en.xm <= en.xp);
  return en;
}

template <class type>
inline eld<type> operator*(const eld<type> e1, const eld<type> e2)
{
  eld<type> en(e1);
  en*= e2;
  
  assert(en.xm <= en.xp);
  return en;
}

template <class type>
inline eld<type> operator/(const eld<type> e1, const eld<type> e2)
{
  eld<type> en(e1);
  en/= e2;
  
  assert(en.xm <= en.xp);
  return en;
}

template <class type>
inline eld<type> operator+(const type l, const eld<type> e)
{
  eld<type> en(l);
  en+= e;

  assert(en.xm <= en.xp);
  return en;
}

template <class type>
inline eld<type> operator-(const type l, const eld<type> e)
{
  eld<type> en(l);
  en-= e;
  
  assert(en.xm <= en.xp);
  return en;
}

template <class type>
inline eld<type> operator*(const type l, const eld<type> e)
{
  eld<type> en(l);
  en*= e;
  
  assert(en.xm <= en.xp);
  return en;
}

template <class type>
inline eld<type> operator/(const type l, const eld<type> e)
{
  eld<type> en(l);
  en/= e;

  assert(en.xm <= en.xp);
  return en;
}

template <class type>
inline eld<type> operator+(const eld<type> e, const type l)
{
  eld<type> en(e);
  en+= l;

  assert(en.xm <= en.xp);
  return en;
}

template <class type>
inline eld<type> operator-(const eld<type> e, const type l)
{
  eld<type> en(e);
  en-= l;

  assert(en.xm <= en.xp);
  return en;
}

template <class type>
inline eld<type> operator*(const eld<type> e, const type l)
{
  eld<type> en(l);
  en*= e;

  assert(en.xm <= en.xp);
  return en;
}

template <class type>
inline eld<type> operator/(const eld<type> e, const type l)
{
  eld<type> en(one/l);
  en*= e;

  assert(en.xm <= en.xp);
  return en;
}


template <class type>
ostream& operator<<(ostream& os, const eld<type> e)
{
  os << "(" << e.xm << " " << e.x << " " << e.xp << ")";
  return os;
}
    

template <class type>
istream& operator>>(istream& is, eld<type>& e)
{
  char buf= ' ';
  while(buf != '(')
  {
    is >> buf;
  }
  is >> e.xm;
  is >> e.x;
  is >> e.xp;
  
  is >> buf;
  
  assert(e.xm <= e.xp);
  return is;
}

template <class type>
inline int operator==(const eld<type> e1, const eld<type> e2)
{
  return e1.x == e2.x;
}

template <class type>
inline int operator!=(const eld<type> e1, const eld<type> e2)
{
  return !(e1 == e2);
}

template <class type>
inline int operator<(const eld<type> e1, const eld<type> e2)
{
  return e1.x < e2.x;
}

template <class type>
inline int operator<=(const eld<type> e1, const eld<type> e2)
{
  return e1.x <= e2.x;
}

template <class type>
inline int operator>(const eld<type> e1, const eld<type> e2)
{
  return e1.x > e2.x;
}

template <class type>
inline int operator>=(const eld<type> e1, const eld<type> e2)
{
  return e1.x >= e2.x;
}


template <class type>
inline int operator==(const eld<type> e, const type l)
{
  return e.x == l;
}

template <class type>
inline int operator!=(const eld<type> e, const type l)
{
  return !(e == l);
}

template <class type>
inline int operator<(const eld<type> e, const type l)
{
  return e.x < l;
}

template <class type>
inline int operator<=(const eld<type> e, const type l)
{
  return e.x <= l;
}

template <class type>
inline int operator>(const eld<type> e, const type l)
{
  return e.x > l;
}

template <class type>
inline int operator>=(const eld<type> e, const type l)
{
  return e.x >= l;
}


template <class type>
inline int operator==(const type l, const eld<type> e)
{
  return e.x == l;
}

template <class type>
inline int operator!=(const type l, const eld<type> e)
{
  return !(e == l);
}

template <class type>
inline int operator<(const type l, const eld<type> e)
{
  return l < e.x;
}

template <class type>
inline int operator<=(const type l, const eld<type> e)
{
  return l <= e.x;
}

template <class type>
inline int operator>(const type l, const eld<type> e)
{
  return l > e.x;
}

template <class type>
inline int operator>=(const type l, const eld<type> e)
{
  return l >= e.x;
}

template <class type>
inline eld<type> min(const eld<type> e, const type l)
{
  if (e < l) return e;
  else
  {
    eld<type> en(l);
    return en;
  }
}

template <class type>
inline eld<type> min(const type l, const eld<type> e)
{
  if (l < e)
  {
    eld<type> en(l);
    return en;
  }
  else return e;
}

template <class type>
inline eld<type> max(const eld<type> e, const type l)
{
  if (e > l) return e;
  else
  {
    eld<type> en(l);
    return en;
  }
}


template <class type>
inline eld<type> max(const type l, const eld<type> e)
{
  if (l > e)
  {
    eld<type> en(l);
    return en;
  }
  else return e;
}


template <class type>
inline eld<type> sqrt(const eld<type> e)
{
  eld<type> en;

  assert(en.x >= null);
  if (e.xm < null) en.xm= null;
  else en.xm= sqrt(e.xm);
  en.x= sqrt(e.x);
  en.xp= sqrt(e.xp);
  en.xm-= abs(en.xm*__eld_eps);
  en.xp+= abs(en.xp*__eld_eps);

  assert(en.xm <= en.xp);
  return en;
}


template <class type>
inline eld<type> sinh(const eld<type> e)
{
  eld<type> en;
  
  en.xm= sinh(e.xm);
  en.x= sinh(e.x);
  en.xp= sinh(e.xp);
  en.xm-= abs(en.xm*__eld_eps);
  en.xp+= abs(en.xp*__eld_eps);
  
  assert(en.xm <= en.xp);
  return en;
}


template <class type>
inline eld<type> cosh(const eld<type> e)
{
  eld<type> en;
  type h1, h2;
  
  h1= cosh(e.xm);
  h2= cosh(e.xp);
  if ((e.xm < 0) && (e.xp > 0)) en.xm= null;
  else en.xm= min(h1, h2);
  
  en.xp= max(h1, h2);
  en.x= cosh(e.x);
  en.xm-= abs(en.xm*__eld_eps);
  en.xp+= abs(en.xp*__eld_eps);
  
  assert(en.xm <= en.xp);  
  return en;
}

template <class type>
inline eld<type> tanh(const eld<type> e)
{
  eld<type> en;
  
  en.xm= tanh(e.xm);
  en.x= tanh(e.x);
  en.xp= tanh(e.xp);
  en.xm-= abs(en.xm*__eld_eps);
  en.xp+= abs(en.xp*__eld_eps);    
  
  assert(en.xm <= en.xp);  
  return en;
}


template <class type>
inline eld<type> log(const eld<type> e)
{
  eld<type> en;
  
  en.xm= log(e.xm);
  en.x= log(e.x);
  en.xp= log(e.xp);
  en.xm-= abs(en.xm*__eld_eps);
  en.xp+= abs(en.xp*__eld_eps);
  
  assert(en.xm <= en.xp);  
  return en;
}


template <class type>
inline eld<type> exp(const eld<type> e)
{
  eld<type> en;

  en.xm= exp(e.xm);
  en.x= exp(e.x);
  en.xp= exp(e.xp);
  en.xm-= abs(en.xm*__eld_eps);
  en.xp+= abs(en.xp*__eld_eps);
  
  assert(en.xm <= en.xp);  
  return en;
}


template <class type>
inline eld<type> atan(const eld<type> e)
{
  eld<type> en;
  
  en.xm= atan(e.xm);
  en.x= atan(e.x);
  en.xp= atan(e.xp);
  en.xm-= abs(en.xm*__eld_eps);
  en.xp+= abs(en.xp*__eld_eps);

  assert(en.xm <= en.xp);  
  return en;
}

/*
template <class type>
inline eld<type> cos(const eld<type> e)
{
    eld<type> en;

    int xmi= (int) floor(e.xm/__eld_pi);
    int xpi= (int) floor(e.xp/__eld_pi);
    if (xmi != xpi)
    {
	if (xpi%2 == 0)
	{
	    en.xm= min(cos(e.xm), cos(e.xp));
	    en.xp= one;
	}
	else
	{
	    en.xm= -one;
	    en.xp= max(cos(e.xm), cos(e.xp));
	}
    }
    else
    {
	en.xm= min(cos(e.xm), cos(e.xp));
	en.xp= max(cos(e.xm), cos(e.xp));
    }
    en.x= cos(e.x);
    en.xm-= abs(en.xm*__eld_eps);
    en.xp+= abs(en.xp*__eld_eps);

    return en;
}


template <class type>
inline eld<type> sin(const eld<type> e)
{
    eld<type> en;

    int xmi= (int) floor((e.xm-__eld_pi/2)/__eld_pi);
    int xpi= (int) floor((e.xp-__eld_pi/2)/__eld_pi);
    if (xmi != xpi)
    {
	if (xpi%2 == 0)
	{
	    en.xm= min(sin(e.xm), sin(e.xp));
	    en.xp= one;
	}
	else
	{
	    en.xm= -one;
	    en.xp= max(sin(e.xm), sin(e.xp));
	}
    }
    else
    {
	en.xm= min(sin(e.xm), sin(e.xp));
	en.xp= max(sin(e.xm), sin(e.xp));
    }
    en.x= sin(e.x);
    en.xm-= abs(en.xm*__eld_eps);
    en.xp+= abs(en.xp*__eld_eps);

    return en;
}
*/


template <class type>
inline eld<type> tan(const eld<type> e)
{
  eld<type> en;
  
  en.xm= tan(e.xm);
  en.x= tan(e.x);
  en.xp= tan(e.xp);
  en.xm-= abs(en.xm*__eld_eps);
  en.xp+= abs(en.xp*__eld_eps);
  
  assert(en.xm <= en.xp);  
  return en;
}

template <class type>
inline eld<type> abs(const eld<type> e)
{
  eld<type> en;

  en.xm= abs(e.xm);
  en.x= abs(e.x);
  en.xp= abs(e.xp);
  if (en.xm > en.xp)
  {
    type zws= en.xm;
    en.xm= en.xp;
    en.xp= zws;
  }

  assert(en.xm <= en.xp);  
  return en;
}

template <class type>
inline eld<type> pow(const eld<type> e, const eld<type> ee)
{
  assert(e > null);

  return exp(ee*log(e));
}

template <class type>
inline eld<type> pow(const eld<type> e, const type ee)
{
  assert(e > null);

  return exp(ee*log(e));
}

template <class type>
inline eld<type> pow(const eld<type> e, const int ee)
{
  assert(e > null);

  eet= (type) ee;
  return exp(eet*log(e));
}

template <class type>
inline eld<type> pow(const type e, const eld<type> ee)
{
  assert(e > null);

  return exp(ee*log(e));
}

template <class type>
inline eld<type> atanh(const eld<type> e)
{
  eld<type> en;

  assert(e > -one);
  assert(e < one);

  en.xm= atanh(e.xm);
  en.x= atanh(e.x);
  en.xp= atanh(e.xp);
  en.xm-= abs(en.xm*__eld_eps);
  en.xp+= abs(en.xp*__eld_eps);
  
  assert(en.xm <= en.xp);  
  return en;
}
    

template <class type>
inline eld<type> floor(const eld<type> e)
{
  eld<type> en;

  en.xm= floor(e.xm);
  en.x= floor(e.x);
  en.xp= floor(e.xp);

  assert(en.xm <= en.xp);  
  return en;
}


template <class type>
inline eld<type> trunc(const eld<type> e)
{
  eld<type> en;

  en.xm= (type)((long)(e.xm));
  en.x= (type)((long)(e.x));
  en.xp= (type)((long)(e.xp));

  assert(en.xm <= en.xp);  
  return en;
}


template <class type>
inline long conv_to_long(const eld<type> e)
{
  return (long) e.x;
}

template <class type>
inline int conv_to_int(const eld<type> e)
{
  return (int) e.x;
}

template <class type>
inline short conv_to_short(const eld<type> e)
{
  return (short) e.x;
}






