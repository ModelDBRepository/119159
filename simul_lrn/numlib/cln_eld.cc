//--------------------------------------------------------------------------
// Author: Thomas Nowotny
//
// Institute: Institut fuer Theoretische Physik
//            Augustusplatz 10-11
//            04109 Leipzig
//
// email to:  nowotny@itp.uni-leipzig.de
//
// initial version: 11/00
// last change: 11/00
//--------------------------------------------------------------------------

#include "cln_eld.h"

template <>
inline eld<cl_F> log(const eld<cl_F> e)
{
    eld<cl_F> en;

    en.xm= As(cl_F)(ln(e.xm));
    en.xm-= abs(en.xm*__eld_eps);
    en.x= As(cl_F)(ln(e.x));
    en.xp= As(cl_F)(ln(e.xp));
    en.xp+= abs(en.xp*__eld_eps);

    return en;
}

template <>
inline eld<cl_F> atanh(const eld<cl_F> e)
{
  eld<cl_F> en;

  assert(e > -one);
  assert(e < one);

  en.xm= As(cl_F)(atanh(e.xm));
  en.xm-= abs(en.xm*__eld_eps);
  en.x= As(cl_F)(atanh(e.x));
  en.xp= As(cl_F)(atanh(e.xp));
  en.xp+= abs(en.xp*__eld_eps);
 
  return en;
}


template <>
inline eld<cl_F> trunc(const eld<cl_F> e)
{
  eld<cl_F> en;

  en.xm= ftruncate(e.xm);
  en.x= ftruncate(e.x);
  en.xp= ftruncate(e.xp);

  return en;
}

template <>
inline eld<cl_F> pow(const eld<cl_F> e, const int ee)
{
  return pow(e, cl_float(cl_I(ee), __eld_eps));
}

template <>
inline long conv_to_long(const eld<cl_F> e)
{
  return cl_I_to_long(truncate1(e.x));
}

template <>
inline int conv_to_int(const eld<cl_F> e)
{
  return cl_I_to_int(truncate1(e.x));
}

template <>
inline short conv_to_short(const eld<cl_F> e)
{
  return (short) cl_I_to_int(truncate1(e.x));
}

