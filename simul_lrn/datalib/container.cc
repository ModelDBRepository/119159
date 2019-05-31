/***********************************************************************/
/* file container.cc                                                   */
/* contains the implementation of class members of class               */
/* container.                                                          */
/***********************************************************************/

#include "container.h"

/***********************************************************************/
/* implementation of class container                                   */
/***********************************************************************/


/***********************************************************************/
/* constructor                                                         */
/***********************************************************************/

template <class type>
container<type>::container()
{
  sz= 0;
}
 
/***********************************************************************/
/* copy constructor                                                    */
/***********************************************************************/

template <class type>
container<type>::container(const container<type>& c)
{
  sz= c.sz;
}

/***********************************************************************/
/* public member operator=                                             */
/* copies the size from the passed container.                          */
/***********************************************************************/

template <class type>
const container<type>& container<type>::operator=(const container<type>& c)
{
  sz= c.size();
  return *this;
}


/***********************************************************************/
/* public member function empty() takes no argument and returns 1 if   */
/* the container has size 0, returns 0 otherwise.                      */
/***********************************************************************/

template <class type>
int container<type>::empty()
{
  return (sz == 0);
}


/***********************************************************************/
/* public member function size() takes no argument and returns the     */ 
/* size of the container.                                              */
/***********************************************************************/

template <class type>
int container<type>::size()
{
  return sz;
}



