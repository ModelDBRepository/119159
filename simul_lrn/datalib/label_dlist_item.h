/**********************************************************************
 file label_dlist_item.h                                                   
 contains the class definition of class label_dlist_item.                  
                                                                     
 warning: this class is only a auxiliary class for class label_dlist.      
 NEVER USE THIS CLASS BY ITSELF.                                     
 exception: use with another data structure, then make sure that     
 label_dlist_item.cc is also included.                                     
**********************************************************************/


#ifndef LABEL_DLIST_ITEM_H
#define LABEL_DLIST_ITEM_H

/**********************************************************************
 class definition of class label_dlist_item. label_dlist_item is the  
 building block for all types of label_dlists (doubly linked lists).       
**********************************************************************/

template <class type>
class label_dlist_item: public label_container_item<type>
{
 public:
   label_dlist_item(int, const type&);                   // constructor
   ~label_dlist_item();                             // destructor
   label_dlist_item<type> *prev;                    // points to predecessor   
   label_dlist_item<type> *next;                    // points to successor
};


#endif
