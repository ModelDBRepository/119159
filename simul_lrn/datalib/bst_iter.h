/***********************************************************************/
/* file bst_iter.h                                                     */
/* contains the class definition of class bst_iter.                    */
/*                                                                     */
/* warning: this class is only a auxiliary class for class base_bst.   */
/* NEVER USE THIS CLASS BY ITSELF.                                     */
/* exception: use with another data structure, then make sure that     */
/* bst_iter.cc is also included.                                       */
/***********************************************************************/


#ifndef BST_ITERATOR_H
#define BST_ITERATOR_H

#include <assert.h>
#include "base_bst.h"
#include "stack.h"

/***********************************************************************/
/* class definition of class bst_iter. bst_iter is an iterator for bst */
/***********************************************************************/

template <class type>
class bst_iterator: public cont_iterator<type>
{
 protected:
   const base_bst<type> *c;                 // bst the iterator works with
   bst_item<type> *cursor;                  // present position in bst c
   stack<bst_item<type> *> stk;             // nodes still to be worked on
                                            // higher up in the tree
 public:
   bst_iterator(const base_bst<type> *);    // constructor
   bst_iterator(const bst_iterator<type>&); // copy constructor
   const bst_iterator<type>& operator=      // copy operator
     (const bst_iterator<type>&);
   virtual ~bst_iterator() { }              // destructor
   void init();                             // initialize
   int step();                              // step to next node
   virtual bst_item<type> *current();       // return the current node
   type c_value();                          // return current value
};

#endif
