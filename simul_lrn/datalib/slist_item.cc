/***********************************************************************/
/* file slist_item.cc                                                  */
/* contains the implementation of class members of class               */
/* slist_item.                                                         */
/***********************************************************************/

#include "slist_item.h"

/***********************************************************************/
/* implementation of class slist_item                                  */
/***********************************************************************/


/***********************************************************************/
/* contructor                                                          */
/***********************************************************************/

template <class type>
slist_item<type>::slist_item(const type& dat):container_item<type>(dat)
{
   next= NULL;
}

/***********************************************************************/
/* copy constructor                                                    */
/***********************************************************************/

template <class type>
slist_item<type>::slist_item(const slist_item<type>& it):
  container_item<type>(it.data)
{
  if (it.next)
  {
    next= new slist_item<type>(*it.next);
  }
  else
  {
    next= NULL;
  }
}

/***********************************************************************/
/* public operator=                                                    */
/* copies all data members of the passed slist_item.                   */
/***********************************************************************/

template <class type>
const slist_item<type>& slist_item<type>::operator=
(const slist_item<type>& it)
{
  this->data= it.data;
  if (next)
  {
    delete next;
  }
  if (it.next)
  {
    next= new slist_item<type>(*it.next);
  }
  else
  {
    next= NULL;
  }
  return *this;
}

/***********************************************************************/
/* destructor                                                          */
/***********************************************************************/

template <class type>
slist_item<type>::~slist_item()
{
}





