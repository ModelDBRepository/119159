/***********************************************************************/
/* file cnt_bst_set.h                                                  */
/* contains the class definition of class cnt_bst_set.                 */
/***********************************************************************/

#ifndef CNT_BST_H
#define CNT_BST_H

template <class type>
class cnt_bst;                                     // forward declaration

#include "base_bst.h"
#include "cnt_bst_item.h"
#include "cnt_bst_iter.h"
#include <limits.h>

/***********************************************************************/
/* class definition of class cnt_bst_set. cnt_bst_set is a set based   */
/* on a binary search tree, which counts multiple elements             */
/***********************************************************************/

template <class type>
class cnt_bst: public base_bst<type>
{
   friend class cnt_bst_iterator<type>;            // bst_iterator may access
                                                   // private data fields
 protected:
   virtual void del_item(bst_item<type> *, 
		 bst_item<type> *);                // workhorse for del
 public:
   cnt_bst() { }                                   // constructor
   cnt_bst(cnt_bst<type>&) { }                     // copy constructor
   const cnt_bst<type>& operator=                  // copy operator
     (const cnt_bst<type>&);
   virtual ~cnt_bst() {};                          // destructor
   virtual cnt_bst_item<type> *add(const type&);   // adds item once
   virtual cnt_bst_item<type> *cnt_add             // adds items
     (const type&, int cnt=1); 
   virtual cnt_bst_item<type> *search(const type&);    // searches item in bst
   virtual int del(const type&, int cnt= INT_MAX); // removes item from bst
   virtual int count(const type&);                 // return count of data
   virtual cnt_bst_iterator<type> *iterator();     // return new iterator
};

#include "cnt_bst_item.cc"
#include "cnt_bst_iter.cc"
#include "cnt_bst.cc"                              // because of templates

#endif
