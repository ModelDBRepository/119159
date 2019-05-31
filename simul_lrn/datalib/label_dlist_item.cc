/**********************************************************************
 file label_dlist_item.cc                                                  
 contains the implementation of class members of class               
 label_dlist_item.                                                         
**********************************************************************/

#include "label_dlist_item.h"

/**********************************************************************
 implementation of class label_dlist_item                                  
**********************************************************************/


/**********************************************************************
 contructor                                                          
**********************************************************************/

template <class type>
label_dlist_item<type>::label_dlist_item(int label, const type& dat):
  label_container_item<type>(label, dat)
{
  prev= NULL;
  next= NULL;
}

/**********************************************************************
 destructor                                                          
**********************************************************************/

template <class type>
label_dlist_item<type>::~label_dlist_item()
{
}

