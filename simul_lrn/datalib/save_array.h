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

#ifndef SARRAY_H
#define SARRAY_H

#include <assert.h>

template <class type>
class sarray
{
 public:
  type *data;
  int size;
  sarray();
  sarray(int);
  ~sarray();
  type& operator[](int);
  void resize(int);
};

#include "save_array.cc"

#endif
