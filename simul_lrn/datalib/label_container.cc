/***********************************************************************/
/* file label_container.cc                                                   */
/* contains the implementation of class members of class               */
/* label_container.                                                          */
/***********************************************************************/

#include "label_container.h"

/***********************************************************************/
/* implementation of class label_container                                   */
/***********************************************************************/


/***********************************************************************/
/* constructor                                                         */
/***********************************************************************/

template <class type>
label_container<type>::label_container()
{
  this->sz= 0;
}
 
/***********************************************************************/
/* copy constructor                                                    */
/***********************************************************************/

template <class type>
label_container<type>::label_container(const label_container<type>& c)
{
  this->sz= c.sz;
}

/***********************************************************************/
/* public member operator=                                             */
/* copies the size from the passed label_container.                          */
/***********************************************************************/

template <class type>
const label_container<type>& label_container<type>::operator=(const label_container<type>& c)
{
  this->sz= c.size();
}


/***********************************************************************/
/* public member function empty() takes no argument and returns 1 if   */
/* the label_container has size 0, returns 0 otherwise.                      */
/***********************************************************************/

template <class type>
int label_container<type>::empty()
{
  return (this->sz == 0);
}


/***********************************************************************/
/* public member function size() takes no argument and returns the     */ 
/* size of the label_container.                                              */
/***********************************************************************/

template <class type>
int label_container<type>::size()
{
  return this->sz;
}



