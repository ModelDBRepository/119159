#include "TNmath.h"


int multinomial(int n, tnvector<int> k)
{
  int x= 1;

  for (int i= n; i > 1; i--)
  {
    x*= i;
  }
  for (int j= 0; j < k.dim(); j++)
  {
    for (int l= k[j]; l > 1; l--)
    {
      x/= l;
    }
  }
  return x;
}

