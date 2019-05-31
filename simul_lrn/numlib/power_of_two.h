//--------------------------------------------------------------------------
// Author: Thomas Nowotny
//
// Institute: Institut fuer Theoretische Physik
//            Augustusplatz 10-11
//            04109 Leipzig
//
// email to:  nowotny@itp.uni-leipzig.de
//
// initial version: 2/00
// last change: 2/00
//--------------------------------------------------------------------------


#ifndef POWER_OF_TWO_H
#define POWER_OF_TWO_H

template <class type>
class power_of_two
{
 public:
  short LBitSize;
  type *tp;

  power_of_two();
  ~power_of_two();
  type operator[](short);
};

template <class type>
class inv_power_of_two
{
 public:
  short LBitSize;
  type *tp;

  inv_power_of_two();
  ~inv_power_of_two();
  type operator[](short);
};

#include "power_of_two.cc"

#endif

  
