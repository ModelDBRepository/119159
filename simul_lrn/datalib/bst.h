/***********************************************************************/
/* file bst.h                                                          */
/* contains the class definition of class bst.                         */
/***********************************************************************/

#ifndef BST_H
#define BST_H

#include "base_bst.h"

/***********************************************************************/
/* class definition of class bst. bst is a binary search tree          */
/***********************************************************************/

template <class type>
class bst: public base_bst<type>
{
   friend class bst_iterator<type>;              // bst_iterator may access
                                                 // private data fields
public:
   bst() { }                                     // constructor
   bst(const bst<type>& tree):                   // copy constructor
     base_bst<type>(tree) { }
   const bst<type>& operator=(const bst<type>&); // copy operator
   virtual ~bst() { }                            // destructor
   virtual bst_item<type> *add(const type&);     // adds item to bst
   virtual bst_item<type> *search(const type&);  // searches item in bst
};

#include "bst.cc"                                // because of templates

#endif
