//-----------------------------------------------------------------
// multiple loop with the restriction that the inner loop vars
// are > the outer ones.
//-----------------------------------------------------------------

#include "multi_loop.h"

template<class type> multi_loop<type>::multi_loop()
{
  min= 0;
  max= 0;
  dpt= 1;
  cnt= new type[dpt];
}

template<class type> multi_loop<type>::multi_loop(type mini, type maxi,
						  short dpti)
{
  min= mini;
  max= maxi;
  dpt= dpti;
  cnt= new type[dpt];
  reset();
}

template<class type> multi_loop<type>::~multi_loop()
{
  delete[] cnt;
}

template<class type> void multi_loop<type>::reset()
{
  cnt[0]= min;
  for(short d= 1; d < dpt; d ++)
  {
    cnt[d]= cnt[d-1]+1;
  }
}

  
template<class type> int multi_loop<type>::step()
{
  d= dpt-1;
  while(1)
  {
    if (cnt[d] < max)
    {
      cnt[d]++;
      if (d < dpt-1)
      {
	d++;
	cnt[d]= cnt[d-1];
      }
      else return 1;
    }
    else
    {
      d--;
      if (d < 0) return 0;
    }
  }
}


template<class type> type multi_loop<type>::operator[](short n)
{
  return cnt[n];
}


template<class type> type *multi_loop<type>::vals()
{
  return cnt;
}



template<class type> short multi_loop<type>::depth()
{
  return dpt;
}


template<class type> void multi_loop<type>::init(type mini, type maxi,
						  short dpti)
{
  min= mini;
  max= maxi;
  dpt= dpti;
  delete[] cnt;
  cnt= new type[dpt];
  reset();
}
