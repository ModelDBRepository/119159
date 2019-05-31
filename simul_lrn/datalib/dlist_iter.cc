/***********************************************************************/
/* file dlist_iter.cc                                                  */
/* contains the implementation of class members of class               */
/* dlist_iterator.                                                     */
/***********************************************************************/

#include "dlist_iter.h"

/***********************************************************************/
/* implementation of class dlist_iterator                              */
/***********************************************************************/


/***********************************************************************/
/* contructor                                                          */
/***********************************************************************/

template <class type>
dlist_iterator<type>::dlist_iterator(const base_dlist<type> *ls):
cont_iterator<type>(ls)
{
  c= ls;
  init();
}

/***********************************************************************/
/* copy constructor                                                    */
/***********************************************************************/

template <class type>
dlist_iterator<type>::dlist_iterator
(const dlist_iterator<type>& sl_iter):cont_iterator<type>(sl_iter)
{
  c= sl_iter.c;
  next= sl_iter.next;
  cursor= sl_iter.cursor;
}


/***********************************************************************/
/* public member operator=                                             */
/* copies all data members from the passed dlist_iterator.             */
/***********************************************************************/

template <class type>
const dlist_iterator<type>& dlist_iterator<type>::operator=
(const dlist_iterator<type>& sl_iter)
{
  c= sl_iter.c;
  next= sl_iter.next;
  cursor= sl_iter.cursor;
}


/***********************************************************************/
/* public member function init                                         */
/* takes no argument, initializes the iterator and sets the pointer    */
/* next to the head and cursor to NULL. thus a first call of member    */
/* function step is necessary before use of current or c_value.        */
/***********************************************************************/

template <class type>
void dlist_iterator<type>::init()
{
   next= c->head;
   cursor= NULL;
}

/***********************************************************************/
/* public member function step                                         */
/* takes no argument and moves cursor and next further in the dlist.   */
/* It returns 1 if the end wasn't reached, 0 otherwise. If the end was */
/* reached, cursor keeps pointing to the last item of the dlist,       */
/* not to NULL.                                                        */
/***********************************************************************/

template <class type>
int dlist_iterator<type>::step()
{
   if (next != NULL)
   {
      cursor= next;
      next=next->next;
      return 1;
   }
   else
   {
      return 0;
   }
}


/***********************************************************************/
/* public member function step_back                                    */
/* takes no argument and moves cursor and next back in the dlist.      */
/* It returns 1 if the head wasn't reached, 0 otherwise. If the head   */
/* was reached, cursor keeps pointing to the first item of the dlist,  */
/* not to NULL.                                                        */
/***********************************************************************/

template <class type>
int dlist_iterator<type>::step_back()
{
  if (cursor->prev != NULL)
  {
    next= cursor;
    cursor= cursor->prev;
    return 1;
  }
  else
  {
    return 0;
  }
}


/***********************************************************************/
/* public member function current                                      */
/* takes no argument and returns a pointer to the current item in the  */
/* dlist (the item cursor points to). If the list was freshly          */
/* initialized, NULL is returned.                                      */
/***********************************************************************/

template <class type>
dlist_item<type> *dlist_iterator<type>::current()
{
   return cursor;
}


/***********************************************************************/
/* public member function c_value                                      */
/* returns the value of the current item in the dlist (the item cursor */
/* points to). Requires cursor to be not NULL.                         */
/***********************************************************************/

template <class type>
type dlist_iterator<type>::c_value()
{
  assert(cursor != NULL);
  return cursor->data;
}

