/***********************************************************************/
/* file slist_iter.h                                                   */
/* contains the class definition of class slist_iterator.              */
/*                                                                     */
/* warning: this class is only a auxiliary class for class base_slist. */
/* NEVER USE THIS CLASS BY ITSELF.                                     */
/* exception: use with another data structure, then make sure that     */
/* slist_iter.cc is also included.                                     */
/***********************************************************************/

#ifndef SLIST_ITERATOR_H
#define SLIST_ITERATOR_H

#include <assert.h>
#include "container_iter.h"

/***********************************************************************/
/* class definition of class slist_iterator. slist_iterator is the     */
/* iterator for slist, stack, slset, etc.                              */
/***********************************************************************/

template <class type>
class slist_iterator: public cont_iterator<type>
{
protected:
   const base_slist<type> *c;                  // slist the iterator works with
   slist_item<type> *cursor;                   // present position in slist c
   slist_item<type> *next;                     // next node after cursor 
                                               // (allows deletion of current
                                               // node)
public:
   slist_iterator(const base_slist<type>*);    // constructor
   slist_iterator(const slist_iterator<type>&);// copy constructor
   const slist_iterator<type>& operator=       // copy operator
     (const slist_iterator<type>&);
   ~slist_iterator() { }                       // destructor (empty)
   void init();                                // initialize
   int step();                                 // step to next node
   slist_item<type> *current();                // return the current node
   type c_value();                             // return current value
   void set_current(slist_item<type> *);            // sets the cursor and next
};

#endif
