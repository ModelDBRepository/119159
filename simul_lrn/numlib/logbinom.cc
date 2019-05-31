/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2002-04-22
  
--------------------------------------------------------------------------*/

long double logbinomial(long double n, long double k)
{
  assert(n >= 0.0L);
  assert(k >= 0.0L);
  assert(k <= n);
  
  long double x= 0.0L;
  long double bi= 0.0L;

  while (k > x)
  {
    bi+= logl(n-x);
    bi-= logl(k-x);
    x++;
  }
  return bi;
}


long double p_binomial_log(long double n, long double p, long double k)
{
  assert(n >= 0.0L);
  assert(k >= 0.0L);
  assert(k <= n);

  if (fabsl(p-1.0) < 1.0e-9L) {
    if (fabsl(k-n) < 1.0e-9L) return 1.0L;
    else return 0.0L;
  }
  
  if (p < 1.0e-9L) {
    if (k < 1.0e-9L) return 1.0L;
    else return 0.0L;
  }

  long double x= 0.0L;
  long double bi= 0.0L;

  while (k > x)
  {
    bi+= logl(n-x);
    bi-= logl(k-x);
    x++;
  }

  bi+= k*logl(p);
  bi+= (n-k)*logl(1.0L-p);
  
  return expl(bi);
}


