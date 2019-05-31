/***********************************************************************/
/* file bst.h                                                          */
/* contains the class definition of class bst.                         */
/***********************************************************************/

#ifndef BASE_BST_H
#define BASE_BST_H

template <class type>
class base_bst;                                // forward declaration

#include <assert.h>
#include "container.h"
#include "bst_item.h"
#include "bst_iter.h"

/***********************************************************************/
/* class definition of class bst. bst is a binary search tree          */
/***********************************************************************/

template <class type>
class base_bst: public container<type>
{
   friend class bst_iterator<type>;              // bst_iterator may access
                                                 // private data fields
protected:
   virtual void del_item(bst_item<type> *, 
		 bst_item<type> *);              // workhorse for del
   bst_item<type> *parent_search(const type&);   // internal search routine
   bst_item<type> *child(bst_item<type> *,       // searches child with
			 const type&);           // the data passed as #2
   void add_item(bst_item<type> *,               // add item #1 at parent #2
		 bst_item<type> *);
   bst_item<type> *head;                         // points to start of bst

public:
   base_bst();                                   // constructor
   base_bst(const base_bst<type>&);
   const base_bst<type>& operator=               // copy operator
     (const base_bst<type>&);
   virtual ~base_bst();                          // destructor
   virtual int in(const type&);                  // is arg in the bst ?
   virtual int del(const type&);                 // removes arg. from bst
   virtual void clear();                         // empties the bst
   virtual bst_iterator<type> *iterator();       // returns a new iterator 
};

#include "bst_item.cc" 
#include "bst_iter.cc" 
#include "base_bst.cc"                           // because of templates

#endif
