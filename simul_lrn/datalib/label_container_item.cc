/***********************************************************************/
/* file label_container_item.cc                                              */
/* contains the implementation of class members of class               */
/* label_container_item.                                                     */
/***********************************************************************/

#include "label_container_item.h"

/***********************************************************************/
/* implementation of class label_container_item                              */
/***********************************************************************/

/***********************************************************************/
/* constructor                                                         */
/***********************************************************************/

template <class type>
label_container_item<type>::label_container_item(int lab, const type& dat)
{
  label= lab;
  data= dat;
}

/***********************************************************************/
/* copy constructor                                                    */
/***********************************************************************/

template <class type>
label_container_item<type>::label_container_item(const label_container_item<type>& it)
{
  label= it.label;
  data= it.data;
}

/***********************************************************************/
/* public member operator=                                             */
/* copies all data members of the passed label_container_item.               */
/***********************************************************************/

template <class type>
const label_container_item<type>& label_container_item<type>::operator=
(const label_container_item<type>& it)
{
  label= it.label;
  data= it.data;
}

