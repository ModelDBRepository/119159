/***********************************************************************/
/* file slist.cc                                                       */
/* contains the implementation of class members of class               */
/* slist.                                                              */
/***********************************************************************/

#include "slist.h"

/***********************************************************************/
/* implementation of class slist                                       */
/***********************************************************************/

/***********************************************************************/
/* operator=                                                           */
/* copies all data members of the passed slist. The resulting slist    */
/* is a new copy of the passed list.                                   */
/***********************************************************************/

template <class type>
const slist<type>& slist<type>::operator=(const slist<type>& sl)
{
  this->sz= sl.sz;
  if (sl.head != NULL)
  {
    this->head= new slist_item<type>(*sl.head);
    this->tail= this->head;
    while (this->tail->next)
    {
      this->tail= this->tail->next;
    }
  }
  else
  {
    if (this->head != NULL)
    {
      delete this->head;
    }
    this->head= NULL;
    this->tail= NULL;
  }
  return *this;
}


/***********************************************************************/
/* public member function prepend                                      */
/* takes one argument of type and prepends an item to the list         */
/* containing the passed data using the more basic member function     */
/* add().                                                              */
/***********************************************************************/

template <class type>
slist_item<type> *slist<type>::prepend(const type& dat)
{
  return add(dat);
}


/***********************************************************************/
/* public member function del_item                                     */
/* takes a pointer on an slist_item and deletes it deletes from the    */
/* slist. Its data is returned.                                        */
/***********************************************************************/

template <class type>
type slist<type>::del_item(slist_item<type>* it)
{
  assert(it != NULL);
  type data= it->data;
  del(data);
  return data;
}


/***********************************************************************/
/* public member function first                                        */
/* takes no argument and returns a pointer to the first item in the    */
/* slist.                                                              */
/***********************************************************************/

template <class type>
slist_item<type> *slist<type>::first()
{
   return this->head;
}

/***********************************************************************/
/* public member function last                                         */
/* takes no argument and returns a pointer to the last item in the     */
/* slist.                                                              */
/***********************************************************************/

template <class type>
slist_item<type> *slist<type>::last()
{
  return this->tail;
}

/***********************************************************************/
/* public memeber function search                                      */
/* takes an argument of type and searches an item containing this data.*/
/* A pointer on the item containing the correct data ist returned if   */
/* found; otherwise NULL is returned.                                  */
/***********************************************************************/

template <class type>
slist_item<type> *slist<type>::search(const type& dat)
{
  if (this->head != NULL)
  {
    slist_item<type> *parent= parent_search(dat);
    if (parent)
    {
      return parent->next;
    }
    else
    {
      return this->head;
    }
  }
  else
  {
    return NULL;
  }
}


/***********************************************************************/
/* public member function succ                                         */
/* takes a pointer to an item in the slist and returns a pointer to    */
/* the successor in the slist.                                         */
/***********************************************************************/

template <class type>
slist_item<type> *slist<type>::succ(slist_item<type> *it)
{
  assert(it != NULL);
  return it->next;
}


/***********************************************************************/
/* public member function insert                                       */
/* takes a pointer to a slist node and an argument type. It inserts a  */
/* new item with the correct data behind the node received as first    */
/* argument.                                                           */ 
/***********************************************************************/

template <class type>
slist_item<type> *slist<type>::insert(slist_item<type> *parent, 
				      const type& dat)
{
  assert(parent != NULL);
  slist_item<type> *it= new slist_item<type>(dat);
  it->next= parent->next;
  parent->next= it;
  this->sz++;
  if (this->tail == parent)
  {
    this->tail= it;
  }
  return it;
}


/***********************************************************************/
/* public member function append                                       */
/* takes an argument of type and appends an item to the slist          */
/* containing the passed data.                                         */
/***********************************************************************/

template <class type>
slist_item<type> *slist<type>::append(const type& dat)
{
  if (this->head != NULL)          // same as if (tail != NULL)         
  {
    return insert(this->tail, dat);
  }
  else
  {
    return add(dat);
  }
}

template <class type>
ostream &operator<<(ostream& os, slist<type>& sl)
{
  slist_iterator<type> *iter= sl.iterator();

  while(iter->step())
  {
    os << iter->c_value();
    os << endl;
  }
  delete iter;

  return os;
}


template <class type>
istream &operator>>(istream& is, slist<type>& sl)
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

