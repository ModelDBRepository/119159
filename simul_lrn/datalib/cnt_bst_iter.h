
/***********************************************************************/
/* file cnt_bst_iter.h                                                 */
/* contains the class definition of class cnt_bst_iter.                */
/*                                                                     */
/* warning: this class is only a auxiliary class for class cnt_bst.    */
/* NEVER USE THIS CLASS BY ITSELF.                                     */
/* exception: use with another data structure, then make sure that     */
/* cnt_bst_iter.cc is also included.                                   */
/***********************************************************************/


#ifndef CNT_BST_ITERATOR_H
#define CNT_BST_ITERATOR_H

#include <assert.h>
#include "cnt_bst.h"

/***********************************************************************/
/* class definition of class bst_iter. bst_ite is an iterator for bst's*/
/***********************************************************************/

template <class type>
class cnt_bst_iterator: public bst_iterator<type>
{
 public:
   cnt_bst_iterator(const cnt_bst<type> *tree):     // constructor
     bst_iterator<type>(tree) { }                 
   cnt_bst_iterator                                 // copy constructor
     (const cnt_bst_iterator<type>& c_i):
       bst_iterator<type>(c_i) { }
   const cnt_bst_iterator<type>& operator=
     (const cnt_bst_iterator<type>& c_i);           // copy operator
   ~cnt_bst_iterator() { }                          // destructor
   cnt_bst_item<type> *current();                   // return current node
   int c_count();                                   // return current count
};

#endif
