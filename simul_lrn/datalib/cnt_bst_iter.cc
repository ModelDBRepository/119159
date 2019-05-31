/***********************************************************************/
/* file cnt_bst_iter.cc                                                */
/* contains the implementation of class members of class               */
/* cnt_bst_iter.                                                       */
/***********************************************************************/

#include "cnt_bst_iter.h"

/***********************************************************************/
/* implementation of class cnt_bst_iter                                */
/***********************************************************************/

/***********************************************************************/
/* public member operator=                                             */
/* copies all data members of the passed cnt_bst_iterator              */
/***********************************************************************/

template <class type>
const cnt_bst_iterator<type>& cnt_bst_iterator<type>::operator=
  (const cnt_bst_iterator<type>& c_iter)
{
  this->c= c_iter.c;
  this->stk= c_iter.stk;
  this->cursor= c_iter.cursor;
  return *this;
}

/***********************************************************************/
/* public member function current                                      */
/* takes no argument and returns a pointer to the current item in the  */
/* cnt_bst.                                                            */
/***********************************************************************/

template <class type>
cnt_bst_item<type> *cnt_bst_iterator<type>::current()
{
   return (cnt_bst_item<type> *) this->cursor;
}

/***********************************************************************/
/* public member function c_count                                      */
/* takes no argument and returns the count of the current item in the  */
/* cnt_bst.                                                            */
/***********************************************************************/

template <class type>
int cnt_bst_iterator<type>::c_count()
{
  assert(this->cursor != NULL);
  cnt_bst_item<type> *curs= (cnt_bst_item<type> *) this->cursor;
  return curs->count;
}

