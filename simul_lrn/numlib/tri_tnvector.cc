/****************************************************************************/
/* Class implementation tri_vector.cc   (triangle vector)                   */
/****************************************************************************/

#include "tri_tnvector.h"

/****************************************************************************/
/* constructor                                                              */
/****************************************************************************/

template<class type>
tri_vector<type>::tri_vector(int t_size)
{
  size=t_size;
  vec= new type[size*(size+1)/2];
}

/****************************************************************************/
/* copy constructor                                                         */
/****************************************************************************/

template<class type>
tri_vector<type>::tri_vector(const tri_vector<type> &tv)
{
  size= tv.size;
  vec= new type[size*(size+1)/2];
  for(i= 0; i < size*(size+1)/2; i++)
  {
    vec[i]= tv.vec[i];
  }
}

/****************************************************************************/
/* public member operator=                                                  */
/* copies all data members from the passed tri_vector.                      */
/****************************************************************************/

template<class type>
const tri_vector<type> &tri_vector<type>::operator=(const tri_vector<type> &tv)
{
  size= tv.size;
  delete vec;
  vec= new type[size*(size+1)/2];
  for(i= 0; i < size*(size+1)/2; i++)
  {
    vec[i]= tv.vec[i];
  }
  return *this;
}

/****************************************************************************/
/* destructor                                                               */
/****************************************************************************/

template <class type>
tri_vector<type>::~tri_vector()
{
  delete vec;
}

/****************************************************************************/
/* public member function get                                               */
/* takes two arguments of type int and returns the value of the cell        */
/* corresponding to these coordinates.                                      */
/****************************************************************************/

template <class type>
type tri_vector<type>::get(int x, int y)
{
  assert(y < size);

  if (x <= y)
  {
    return vec[y*(y+1)/2+x];
  }
  else
  {
    return -vec[x*(x+1)/2+y];
  }
}

/****************************************************************************/
/* public member function set                                               */
/* takes two coordinates and a value of type and sets the value of the      */
/* cell corresponding to the coordinates to the value passed.               */
/****************************************************************************/

template <class type>
void tri_vector<type>::set(int x, int y, type value)
{
  assert(y < size);

  if (x <= y)
  {
    vec[y*(y+1)/2+x]= value;
  }
  else
  {
    vec[x*(x+1)/2+y]= -value;
  }
}




