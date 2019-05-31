  

template <class Type>
dpolynom<Type>::dpolynom()
{
  dimen= 1;
  ord= 0;
}
 
template <class Type>
dpolynom<Type>::dpolynom(int dim, int order)
{
  dimen= dim;
  ord= order;
  coeff.resize(coeff_size());
  org.resize(dim);
}

template <class Type>
dpolynom<Type>::dpolynom(int dim, int order, tnvector<Type>& a, 
			     tnvector<Type>& origin)
{
  dimen= dim;
  ord= order;
  coeff= a;
  org= origin;
}


template <class Type>
dpolynom<Type>::dpolynom(dpolynom<Type>& p)
{
  dimen= p.dimen;
  ord= p.ord;
  coeff= p.coeff;
  org= p.org;
}

template <class Type>
dpolynom<Type>::dpolynom(istream& is)
{
  read(is);
}

template <class Type>
dpolynom<Type>::~dpolynom()
{
}


template <class Type>
void dpolynom<Type>::resize(int dim, int order)
{
  dimen= dim;
  ord= order;
  coeff.resize(coeff_size());
  org.resize(dim);
}  


template <class Type>
int dpolynom<Type>::dim()
{
  return dimen;
}

template <class Type>
int dpolynom<Type>::order()
{
  return ord;
}

template <class Type>
int dpolynom<Type>::coeff_size()
{
  int sz= 0;
  dclock clk(-dimen, 1, dimen-1, 0, ord-1);
  while(clk.advance())
  {
    sz++;
  }
  clk.reset(-dimen, 1, dimen-1, 0, ord);
  while(clk.advance())
  {
    sz++;
  }
  return sz;
}

template <class Type>
Type dpolynom<Type>::value(tnvector<Type> x)
{
  assert(dimen == x.dim());

  dclock clk;
  int pos= 0, ibuf;
  Type buf;
  Type result= 0;
  x-= org;
  for(int order= ord-1; order <= ord; order++)
  {
    clk.reset(-dimen, 1, dimen-1, 0, order);
    while(clk.advance())
    {
      buf= coeff[pos];
      for(int j= 0; j < order; j++)
      {
	ibuf= clk[j];
	if (ibuf < 0)
	{
	  ibuf= abs(ibuf)-1;
	  if (fabs(x[ibuf]) > epsilon)
	  { 
	    buf/= x[ibuf];
	  }
	  else
	  {
	    return 0;
	  }
	}
	else
	{
	  buf*= x[ibuf];
	}
      }
      result+= buf;
      pos++;
    }
  }
  return result;
}


template <class Type>
Type dpolynom<Type>::coefficient(tnvector<int> i)
{
  assert((i.dim() == ord-1) || (i.dim() == ord));

  int pos= 0;
  dclock clk;
  for (int order= ord-1; order <= ord; order++)
  {
    clk.reset(-dimen, 1, dimen-1, 0, order);
    while(clk.advance())
    {
      if (i == clk.all())
      {
	return coeff[pos];
      }
      else
      {
	pos++;
      }
    }
  }
  return 0;
}

template <class Type>
void dpolynom<Type>::set_coeff_vec(tnvector<Type> co)
{
  assert(coeff_size() == co.dim());

  coeff= co;
}

template <class Type>
void dpolynom<Type>::set_coeff_direct(int pos, Type value)
{
  assert(coeff_size() > pos);

  coeff.set(pos, value);
}


template <class Type>
tnvector<Type> dpolynom<Type>::coeff_vec()
{
  return coeff;
}

template <class Type>
void dpolynom<Type>::set_origin(tnvector<Type> orig)
{
  assert(orig.dim() == dimen);

  org= orig;
}

template <class Type>
tnvector<Type> dpolynom<Type>::origin()
{
  return org;
}

template <class Type>
void dpolynom<Type>::read(istream& is)
{
  char pk;
  pk= is.peek();
  while((pk == '#') || (pk == '\n'))
  {
    is.ignore(80, '\n');
    pk= is.peek();
  }  
  is.ignore(80,' ');
  is >> dimen;
  is.ignore(80,' ');
  is >> ord;
  is.ignore(80,' ');
  org.resize(dimen);
  is >> org;
  pk= is.peek();
  is.ignore(80,'\n');
  coeff.resize(coeff_size());  
  is >> coeff;
}


template <class Type>
void dpolynom<Type>::nice_output(ostream& os)
{
  os << "dimension: " << dimen << endl;
  os << "order: " << ord << endl;
  os << "origin: ";
  os << org;
  os << endl;
  os << "coefficients: " << endl;
  dclock clk;
  int pos= 0;
  for (int order= ord-1; order <= ord; order++)
  {
    clk.reset(-dimen, 1, dimen-1, 0, order);
    while(clk.advance())
    {
      os << clk.all();
      os << ": ";
      os << coeff[pos];
      os << endl;
      pos++;
    }
    os << endl;
  }
  os << endl;
}


template <class Type>
istream& operator>>(istream& is, dpolynom<Type>& p)
{
  p.read(is);
}


template <class Type>
ostream& operator<<(ostream& os, dpolynom<Type>& p)
{
  os << "dimension: " << p.dimen << endl;
  os << "order: " << p.ord << endl;
  os << "origin: ";
  os << p.org;
  os << endl;
  os << "coefficients: " << endl;
  os << p.coeff;
  os << endl;
  return os;
}

