/***********************************************************************/
/* file slist_iter.cc                                                  */
/* contains the implementation of class members of class               */
/* slist_iterator.                                                     */
/***********************************************************************/

#include "slist_iter.h"

/***********************************************************************/
/* implementation of class slist_iterator                              */
/***********************************************************************/


/***********************************************************************/
/* contructor                                                          */
/***********************************************************************/

template <class type>
slist_iterator<type>::slist_iterator(const base_slist<type> *ls):
cont_iterator<type>(ls)
{
  c= ls;
  init();
}

/***********************************************************************/
/* copy constructor                                                    */
/***********************************************************************/

template <class type>
slist_iterator<type>::slist_iterator
(const slist_iterator<type>& sl_iter):cont_iterator<type>(sl_iter)
{
  c= sl_iter.c;
  next= sl_iter.next;
  cursor= sl_iter.cursor;
}


/***********************************************************************/
/* public member operator=                                             */
/* copies all data members from the passed slist_iterator.             */
/***********************************************************************/

template <class type>
const slist_iterator<type>& slist_iterator<type>::operator=
(const slist_iterator<type>& sl_iter)
{
  c= sl_iter.c;
  next= sl_iter.next;
  cursor= sl_iter.cursor;
  return *this;
}


/***********************************************************************/
/* public member function init                                         */
/* takes no argument, initializes the iterator and sets the pointer    */
/* next to the head and cursor to NULL. thus a first call of member    */
/* function step is necessary before use of current or c_value.        */
/***********************************************************************/

template <class type>
void slist_iterator<type>::init()
{
   next= c->head;
   cursor= NULL;
}

/***********************************************************************/
/* public member function step                                         */
/* takes no argument and moves cursor and next further in the slist.   */
/* It returns 1 if the end wasn't reached, 0 otherwise. If the end was */
/* reached, cursor keeps pointing to the last item of the slist,       */
/* not to NULL.                                                        */
/***********************************************************************/

template <class type>
int slist_iterator<type>::step()
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
/* public member function current                                      */
/* takes no argument and returns a pointer to the current item in the  */
/* slist (the item cursor points to). If the list was freshly          */
/* initialized, NULL is returned.                                      */
/***********************************************************************/

template <class type>
slist_item<type> *slist_iterator<type>::current()
{
   return cursor;
}


/***********************************************************************/
/* public member function c_value                                      */
/* returns the value of the current item in the slist (the item cursor */
/* points to). Requires cursor to be not NULL.                         */
/***********************************************************************/

template <class type>
type slist_iterator<type>::c_value()
{
  assert(cursor != NULL);
  return cursor->data;
}

/***********************************************************************/
/* public member function set_current                                  */
/* sets the value of the current an next item in the slist (the item   */
/* cursor and next point to points to). Requires arg to be not NULL    */
/* and better belong to the list.            .                         */
/***********************************************************************/

template <class type>
void slist_iterator<type>::set_current(slist_item<type> *nc)
{
  assert(nc != NULL);
  cursor= nc;
  next= nc->next;
}
  
