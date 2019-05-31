//--------------------------------------------------------------------------
// Author: Thomas Nowotny
//
// Institute: Institut fuer Theoretische Physik
//            Augustusplatz 10-11
//            04109 Leipzig
//
// email to:  nowotny@itp.uni-leipzig.de
//
// initial version: 11/98
// last change: 11/98
//--------------------------------------------------------------------------

#include "save_array.h"


template <class type>
sarray<type>::sarray()
{
  size= 0;
  data= new type[1];
}


template <class type>
sarray<type>::sarray(int sz)
{
  size= sz;
  data= new type[sz];
}

template <class type>
sarray<type>::~sarray()
{
  delete[] data;
}

template <class type>
type& sarray<type>::operator[](int i)
{
  assert(i >= 0);
  assert(i < size);

  return data[i];
}


template <class type>
void sarray<type>::resize(int sz)
{
  assert(sz > 0);

  delete[] data;
  data= new type[sz];
  size= sz;
}

