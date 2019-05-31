/***********************************************************************/
/* file dlist_iter.h                                                   */
/* contains the class definition of class dlist_iterator.              */
/*                                                                     */
/* warning: this class is only a auxiliary class for class base_dlist. */
/* NEVER USE THIS CLASS BY ITSELF.                                     */
/* exception: use with another data structure, then make sure that     */
/* dlist_iter.cc is also included.                                     */
/***********************************************************************/

#ifndef DLIST_ITERATOR_H
#define DLIST_ITERATOR_H

#include <assert.h>

/***********************************************************************/
/* class definition of class dlist_iterator. dlist_iterator is the     */
/* iterator for dlist, stack, slset, etc.                              */
/***********************************************************************/

template <class type>
class dlist_iterator: public cont_iterator<type>
{
protected:
   const base_dlist<type> *c;                  // dlist the iterator works with
   dlist_item<type> *cursor;                   // present position in dlist c
   dlist_item<type> *next;                     // next node after cursor 
                                               // (allows deletion of current
                                               // node)
public:
   dlist_iterator(const base_dlist<type>*);    // constructor
   dlist_iterator(const dlist_iterator<type>&);// copy constructor
   const dlist_iterator<type>& operator=       // copy operator
     (const dlist_iterator<type>&);
   ~dlist_iterator() { }                       // destructor (empty)
   void init();                                // initialize
   int step();                                 // step to next node
   int step_back();                            // step to previous node
   dlist_item<type> *current();                // return the current node
   type c_value();                             // return current value
};

#endif
