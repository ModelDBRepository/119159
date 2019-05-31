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

#include "power_of_two.h"

template <class type>
power_of_two<type>::power_of_two()
{
  LBitSize= sizeof(type)*8;
  tp= new type[LBitSize]; 
  tp[0]= 1;                         
  for (int i= 0; i < LBitSize-1; i++)
  {
    tp[i+1]= tp[i]*2;
  }
}

template <class type>
power_of_two<type>::~power_of_two()
{
  delete[] tp;
}

template <class type>
inline type power_of_two<type>::operator[](short i)
{
  return tp[i];
}

template <class type>
inv_power_of_two<type>::inv_power_of_two()
{
  LBitSize= sizeof(type)*8;
  tp= new type[LBitSize]; 
  tp[LBitSize-1]= 1;                         
  for (int i= LBitSize-1; i > 0; i--)
  {
    tp[i-1]= tp[i]*2;
  }
}

template <class type>
inv_power_of_two<type>::~inv_power_of_two()
{
  delete[] tp;
}

template <class type>
inline type inv_power_of_two<type>::operator[](short i)
{
  return tp[i];
}

power_of_two<ulong> twopow;
inv_power_of_two<ulong> invtp;
power_of_two<int> pow2int;
power_of_two<long> pow2long;
power_of_two<long double> pow2ld;

