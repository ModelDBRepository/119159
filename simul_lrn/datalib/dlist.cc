/***********************************************************************/
/* file dlist.cc                                                       */
/* contains the implementation of class members of class               */
/* dlist.                                                              */
/***********************************************************************/

#include "dlist.h"

/***********************************************************************/
/* implementation of class dlist                                       */
/***********************************************************************/

/***********************************************************************/
/* public member function prepend                                      */
/* takes one argument of type and prepends an item to the list         */
/* containing the passed data using the more basic member function     */
/* add().                                                              */
/***********************************************************************/

template <class type>
dlist_item<type> *dlist<type>::prepend(const type& dat)
{
  return add(dat);
}

/***********************************************************************/
/* public member function first                                        */
/* takes no argument and returns a pointer to the first item in the    */
/* dlist.                                                              */
/***********************************************************************/

template <class type>
dlist_item<type> *dlist<type>::first()
{
   return this->head;
}

/***********************************************************************/
/* public member function last                                         */
/* takes no argument and returns a pointer to the last item in the     */
/* dlist.                                                              */
/***********************************************************************/

template <class type>
dlist_item<type> *dlist<type>::last()
{
  return this->tail;
}

/***********************************************************************/
/* public member function succ                                         */
/* takes a pointer to an item in the dlist and returns a pointer to    */
/* the successor in the dlist.                                         */
/***********************************************************************/

template <class type>
dlist_item<type> *dlist<type>::succ(dlist_item<type> *it)
{
  assert(it != NULL);
  return it->next;
}

template <class type>
dlist_item<type> *dlist<type>::search(const type& dat)
{
  return find(dat);
}

/***********************************************************************/
/* public member function insert                                       */
/* takes a pointer to a dlist node and an argument type. It inserts a  */
/* new item with the correct data behind the node received as first    */
/* argument.                                                           */ 
/***********************************************************************/

template <class type>
dlist_item<type> *dlist<type>::insert_after(dlist_item<type> *parent, 
					    const type& dat)
{
  assert(parent != NULL);
  dlist_item<type> *it= new dlist_item<type>(dat);
  it->next= parent->next;
  it->prev= parent;
  parent->next= it;
  this->sz++;
  if (this->tail == parent) {
    this->tail= it;
  }
  else {
    it->next->prev= it;
  }
  return it;
}


/***********************************************************************/
/* public member function append                                       */
/* takes an argument of type and appends an item to the dlist          */
/* containing the passed data.                                         */
/***********************************************************************/

template <class type>
dlist_item<type> *dlist<type>::append(const type& dat)
{
  if (this->head != NULL)          // same as if (tail != NULL)         
  {
    return insert_after(this->tail, dat);
  }
  else
  {
    return add(dat);
  }
}

template <class type>
ostream &operator<<(ostream& os, dlist<type>& sl)
{
  dlist_iterator<type> *iter= sl.iterator();

  forall (iter)
  {
    os << iter->c_value();
    os << endl;
  }
  delete iter;

  return os;
}


template <class type>
istream &operator>>(istream& is, dlist<type>& sl)
{
  type dat;
  while (is.good())
  {
    is >> dat;
    if (is.good())
    {
      sl.append(dat);
    }
  }
  return is;
}

