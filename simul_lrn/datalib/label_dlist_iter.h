/*********************************************************************
 file label_dlist_iter.h                                                   
 contains the class definition of class label_dlist_iterator.              
                                                                     
 warning: this class is only a auxiliary class for class base_label_dlist. 
 NEVER USE THIS CLASS BY ITSELF.                                     
 exception: use with another data structure, then make sure that     
 label_dlist_iter.cc is also included.                                     
*********************************************************************/

#ifndef LABEL_DLIST_ITERATOR_H
#define LABEL_DLIST_ITERATOR_H

#include <assert.h>

/*********************************************************************
 class definition of class label_dlist_iterator. label_dlist_iterator
 is the iterator for label_dlist, stack, slset, etc.
*********************************************************************/

template <class type>
class label_dlist_iterator: public label_cont_iterator<type>
{
protected:
   const base_label_dlist<type> *c;     // label_dlist the iterator works with
   label_dlist_item<type> *cursor;      // present position in label_dlist c
   label_dlist_item<type> *next;        // next node after cursor 
                                        // (allows deletion of current node)
public:
   label_dlist_iterator(const base_label_dlist<type>*);    // constructor
   label_dlist_iterator(const label_dlist_iterator<type>&);// copy constructor
   const label_dlist_iterator<type>& operator=             // copy operator
     (const label_dlist_iterator<type>&);
   ~label_dlist_iterator() { }                            // destructor (empty)
   void init();                                // initialize
   int step();                                 // step to next node
   label_dlist_item<type> *current();          // return the current node
   type c_value();                             // return current value
   int c_label();
};

#endif
