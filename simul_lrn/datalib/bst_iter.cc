/***********************************************************************/
/* file bst_iter.cc                                                    */
/* contains the implementation of class members of class               */
/* bst_iter.                                                           */
/***********************************************************************/

#include "bst_iter.h"

/***********************************************************************/
/* implementation of class bst_iter                                    */
/***********************************************************************/

/***********************************************************************/
/* constructor                                                         */
/***********************************************************************/

template <class type>
bst_iterator<type>::bst_iterator(const base_bst<type> *tree):
cont_iterator<type>(tree)
{
  assert(tree != NULL);
  c= tree;
  init();
} 

/***********************************************************************/
/* copy constructor                                                    */
/* only attaches the iterator to the same tree but initializes it to   */
/* to the first position in the tree.                                  */ 
/***********************************************************************/

template <class type>
bst_iterator<type>::bst_iterator(const bst_iterator<type>& bst_iter):
cont_iterator<type>(bst_iter)
{
  c= bst_iter.c;
  stk= bst_iter.stk;
  cursor= bst_iter.cursor;
}

/***********************************************************************/
/* public member operator=                                             */
/* copies all values from the passed bst_iterator including the        */
/* complete stack and current position.                                */
/***********************************************************************/

template <class type>
const bst_iterator<type>& bst_iterator<type>::operator=
(const bst_iterator<type>& bst_iter)
{
  c= bst_iter.c;
  stk= bst_iter.stk;
  cursor= bst_iter.cursor;
  return *this;
}

/***********************************************************************/
/* public member function init                                         */
/* takes no arguments and initializes the iterator. it sets the cursor */
/* to the leftmost node in the tree and saves the parent nodes on the  */
/* stack stk. Cursor points beyond the leftmost node to the left (i.e. */
/* is NULL) such that one step() is necessary prior to use of          */
/* member function current() or c_value().                             */
/***********************************************************************/

template <class type>
void bst_iterator<type>::init()
{
  stk.clear();
  cursor= c->head;
  while(cursor)
  {
    stk.push(cursor);
    cursor= cursor->left;
  }
}

/***********************************************************************/
/* public member function step                                         */
/* moves one node further in an inorder traversal of the tree c.       */
/* It returns 0 if the traversal was already completed, 1 otherwise.   */
/***********************************************************************/

template <class type>
int bst_iterator<type>::step()
{
   if (cursor)
   {
     if (cursor->right)
     {
       cursor= cursor->right;
       while(cursor)
       {
	 stk.push(cursor);
	 cursor= cursor->left;
       }
     }
   }
   if (stk.empty())
   {
     return 0;
   }
   else
   {
     cursor= stk.pop();
     return 1;
   }
}


/***********************************************************************/
/* public member function current                                      */
/* returns a pointer to the current item in the bst, i.e. the item     */
/* the cursor currently points to.                                     */
/***********************************************************************/

template <class type>
bst_item<type> *bst_iterator<type>::current()
{
   return cursor;
}

/***********************************************************************/
/* public member function c_value                                      */
/* returns the value of the current item in the bst, i.e. the value of */
/* the item the cursor currently points to. If cursor points to NULL,  */
/* 0 is returned (type needs to have some 0 element).                  */
/***********************************************************************/

template <class type>
type bst_iterator<type>::c_value()
{
  assert(cursor != NULL);
  return cursor->data;
}


