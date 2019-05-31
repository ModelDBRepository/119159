/***********************************************************************/
/* file container_item.cc                                              */
/* contains the implementation of class members of class               */
/* container_item.                                                     */
/***********************************************************************/

#include "container_item.h"

/***********************************************************************/
/* implementation of class container_item                              */
/***********************************************************************/

/***********************************************************************/
/* constructor                                                         */
/***********************************************************************/

template <class type>
container_item<type>::container_item(const type& dat)
{
   data= dat;
}

/***********************************************************************/
/* copy constructor                                                    */
/***********************************************************************/

template <class type>
container_item<type>::container_item(const container_item<type>& it)
{
   data= it.data;
}

/***********************************************************************/
/* public member operator=                                             */
/* copies all data members of the passed container_item.               */
/***********************************************************************/

template <class type>
const container_item<type>& container_item<type>::operator=
(const container_item<type>& it)
{
   data= it.data;
   return *this;
}

