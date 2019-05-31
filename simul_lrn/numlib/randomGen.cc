/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2002-02-01
  
--------------------------------------------------------------------------*/

#include "randomGen.h"

randomGen::randomGen()
{
  time_t t0= time(NULL);
  assert (t0 != -1);

  ulong seed1= t0;
  ulong seed2= t0%131313;
  ulong seed3= t0%171717;
  TheGen.srand(seed1, seed2, seed3);
}

double randomGen::n()
{
  a= TheGen.rand();
  a/= ULONG_MAX;
  return a;
}
  

stdRG::stdRG()
{
  time_t t0= time(NULL);
  srand(t0);
  themax= RAND_MAX+1.0;
}


double stdRG::n()
{
  return rand()/themax;
}
