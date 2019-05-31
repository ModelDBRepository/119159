/*********************************************************************
 file label_dlist.cc                                                       
 contains the implementation of class members of class               
 label_dlist.                                                              
*********************************************************************/

#include "label_dlist.h"

/*********************************************************************
 implementation of class label_dlist                                       
*********************************************************************/

/*********************************************************************
 public member function first                                        
 takes no argument and returns a pointer to the first item in the    
 label_dlist.                                                              
*********************************************************************/

template <class type>
label_dlist_item<type> *label_dlist<type>::first()
{
   return head;
}

/*********************************************************************
 public member function last                                         
 takes no argument and returns a pointer to the last item in the     
 label_dlist.                                                              
*********************************************************************/

template <class type>
label_dlist_item<type> *label_dlist<type>::last()
{
  return tail;
}


/*********************************************************************
 public member function find_left                                        

! not to be used on an empty list!
*********************************************************************/

template <class type>
int label_dlist<type>::find_left(int lb, label_dlist_item<type> *&it)
{
  assert (head != NULL);
  it= head;
  while ((it->next != NULL) && (it->next->label <= lb)) {
    it= it->next;
  }
  return (it->label == lb);
}
