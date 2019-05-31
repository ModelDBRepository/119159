/*********************************************************************
 file label_dlist.h                                                        
 contains the class definition of class label_dlist                        
*********************************************************************/

#ifndef LABEL_DLIST_H
#define LABEL_DLIST_H

#include <assert.h>
#include "base_label_dlist.h"

/*********************************************************************
 class definition of class label_dlist. label_dlist is a singly linked
 list.     
*********************************************************************/

template <class type>
class label_dlist: public base_label_dlist<type>
{
   friend class label_dlist_iterator<type>;  // label_dlist_iterator may access
                                             // private data fields

public:
   label_dlist() { }                                     // constructor
   virtual ~label_dlist() { }                            // destructor (empty)
   virtual label_dlist_item<type> *first();           // returns first node
   virtual label_dlist_item<type> *last();            // returns last node
   virtual int find_left(int, label_dlist_item<type> *&);
};

#include "label_dlist.cc"                            // because of templates

#endif
