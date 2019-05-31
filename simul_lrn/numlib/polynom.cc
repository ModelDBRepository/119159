  

template <class Type>
polynomial<Type>::polynomial()
{
  dimen= 1;
  ord= 0;
}
 
template <class Type>
polynomial<Type>::polynomial(int dim, int order)
{
  dimen= dim;
  ord= order;
  coeff.resize(coeff_size());
  org.resize(dim);
}

template <class Type>
polynomial<Type>::polynomial(int dim, int order, tnvector<Type>& a, 
			     tnvector<Type>& origin)
{
  dimen= dim;
  ord= order;
  coeff= a;
  org= origin;
}


template <class Type>
polynomial<Type>::polynomial(polynomial<Type>& p)
{
  dimen= p.dimen;
  ord= p.ord;
  coeff= p.coeff;
  org= p.org;
}

template <class Type>
polynomial<Type>::polynomial(istream& is)
{
  read(is);
}

template <class Type>
polynomial<Type>::~polynomial()
{
}


template <class Type>
void polynomial<Type>::resize(int dim, int order)
{
  dimen= dim;
  ord= order;
  coeff.resize(coeff_size());
  org.resize(dim);
}  


template <class Type>
int polynomial<Type>::dim()
{
  return dimen;
}

template <class Type>
int polynomial<Type>::order()
{
  return ord;
}

template <class Type>
int polynomial<Type>::coeff_size()
{
  int sz;
  if (dimen > 1)
  {
    sz= (int) pow(dimen, ord+1)-1;
    sz/= dimen-1;
  }
  else
  {
    sz= ord+1;
  }
  return sz;
}

template <class Type>
Type polynomial<Type>::value(tnvector<Type> x)
{
  assert(dimen == x.dim());

  tnvector<int> idx;
  int level;
  Type result= coeff[0];
  Type buf;

  x-= org;
  for(int order= 1; order <= ord; order++)
  {
    idx.resize(order);
    level= 0;
    idx.set(0, -1);
    while(level != -1)
    {
      idx.set(level, idx[level]+1);
      if (idx[level] == dimen)
      {
	level--;
      }
      else
      {
	if (level == order-1)
	{
	  buf= coefficient(idx);
	  for(int j= 0; j < order; j++)
	  {
	    buf*= x[idx[j]];
	  }
	  result+= buf;
	}
	else
	{
	  level++;
	  idx.set(level, -1);
	}
      }
    }
  }
  return result;
}


template <class Type>
Type polynomial<Type>::coefficient(tnvector<int> i)
{
  int pos;
  if (dimen > 1)
  {
    pos= (int) pow(dimen, i.dim())-1;
    pos/= dimen-1;
  }
  else
  {
    pos= i.dim();
  }
  int base= 1;
  for(int k= 0; k < i.dim(); k++)
  {
    pos+= i[i.dim()-k-1]*base;
    base*= dimen;
  }
  return coeff[pos];
}

template <class Type>
void polynomial<Type>::set_coeff_vec(tnvector<Type> co)
{
  assert(coeff_size() == co.dim());

  coeff= co;
}

template <class Type>
void polynomial<Type>::set_coeff_direct(int pos, Type value)
{
  assert(coeff_size() > pos);

  coeff.set(pos, value);
}


template <class Type>
tnvector<Type> polynomial<Type>::coeff_vec()
{
  return coeff;
}

template <class Type>
void polynomial<Type>::set_origin(tnvector<Type> orig)
{
  assert(orig.dim() == dimen);

  org= orig;
}

template <class Type>
tnvector<Type> polynomial<Type>::origin()
{
  return org;
}

template <class Type>
void polynomial<Type>::read(istream& is)
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
istream& operator>>(istream& is, polynomial<Type>& p)
{
  p.read(is);
}

template <class Type>
void polynomial<Type>::nice_output(ostream& os)
{
  os << "dimension: " << dimen << endl;
  os << "order: " << ord << endl;
  os << "origin: ";
  os << org;
  os << endl;
  os << "coefficients: " << endl;
  clock clk;
  for (int order= 0; order <= ord; order++)
  {
    clk.reset(0, 0, dimen-1, 0, order);
    while(clk.advance())
    {
      os << clk.all();
      os << ": ";
      os << coefficient(clk.all());
      os << endl;
    }
    os << endl;
  }
  os << endl;
}

template <class Type>
ostream& operator<<(ostream& os, polynomial<Type>& p)
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
