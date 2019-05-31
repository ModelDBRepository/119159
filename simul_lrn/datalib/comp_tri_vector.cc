/****************************************************************************/
/* Class implementation comp_tri_vector.cc   (triangle vector)              */
/****************************************************************************/

#include "comp_tri_vector.h"

comp_tri_vector::comp_tri_vector(int t_size)
{
  size= t_size;
  int vec_size= size*(size+1)/2;
  vec_size/=4;
  vec_size++;
  vec= new unsigned char[vec_size];
  for(int i= 0; i < vec_size; i++)
  {
    vec[i]= 0;
  }
}


comp_tri_vector::~comp_tri_vector()
{
  delete vec;
}


int comp_tri_vector::get(int x, int y)
{
  int fac= 0;
  int pos, i_pos;
  if (x <= y)
  {
    pos= y*(y+1)/2+x;
    fac= 1;
  }
  else
  {
    pos= x*(x+1)/2+y;
    fac= -1;
  }
  i_pos= pos % 4;
  pos/= 4;
  for (int i= 0; i<=2; i+=2)
  {
    if (vec[pos] & pat[i][i_pos])
    {
      return (i-1)*fac;
    }
  }
  return 0;
}


void comp_tri_vector::set(int x, int y, int value)
{
  int fac= 0;
  int pos, i_pos;
  if (x <= y)
  {
    pos= y*(y+1)/2+x;
    fac= 1;
  }
  else
  {
    pos= x*(x+1)/2+y;
    fac= -1;
  }
  i_pos= pos % 4;
  pos/= 4;
  vec[pos]= ((vec[pos] & pat[3][i_pos]) | pat[value+1][i_pos]) * fac;
}


