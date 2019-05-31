/***********************************************************************/
/* file bst_item.h                                                     */
/* contains the class definition of class bst_item.                    */
/*                                                                     */
/* warning: this class is only a auxiliary class for class bst.        */
/* NEVER USE THIS CLASS BY ITSELF.                                     */
/* exception: use with another data structure, then make sure that     */
/* bst_item.cc is also included.                                       */
/***********************************************************************/

#ifndef BST_ITEM_H
#define BST_ITEM_H

#include "container.h"

/***********************************************************************/
/* class definition of class bst_item. bst_item is the building block  */
/* for binary search trees (bst).                                      */
/***********************************************************************/

template <class type>
class bst_item: public container_item<type>
{
 public:
   bst_item(const type&);                      // constructor
   bst_item(const bst_item<type>&);            // copy constructor
   const bst_item<type>& operator=             // copy operator
     (const bst_item<type>&);
   virtual ~bst_item();                        // destructor
   bst_item<type> *left;                       // point to successors
   bst_item<type> *right;
};

#endif
