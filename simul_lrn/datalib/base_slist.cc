/***********************************************************************/
/* file base_slist.cc                                                  */
/* contains the implementation of class members of class               */
/* base_slist.                                                         */
/***********************************************************************/

#include "base_slist.h"


/***********************************************************************/
/* implementation of class base_slist                                  */
/***********************************************************************/

/***********************************************************************/
/* protected member function internal_search()                         */
/* takes an argument of type and searches for the first item in the    */
/* slist with identical data as the passed argument. Then it returns   */
/* a pointer to the the item in front of the item with the relevant    */
/* data. If the data is not found, a pointer on the last node is       */
/* returned                                                            */
/***********************************************************************/

template <class type>
slist_item<type> *base_slist<type>::parent_search(const type& dat)
{
  slist_item<type> *parent= NULL;
  slist_item<type> *it= head;

  while (it)
  {
    if (it->data == dat)
    {
      return parent;
    }
    parent= it;
    it= it->next;
  }
  return parent;
}


/***********************************************************************/
/* constructor                                                         */
/***********************************************************************/

template <class type>
base_slist<type>::base_slist()
{
   head = NULL;
   tail= NULL;
}

/***********************************************************************/
/* copy constructor                                                    */
/***********************************************************************/

template <class type>
base_slist<type>::base_slist(const base_slist<type>& sl)
{
  this->sz= sl.size();
  if (sl.head != NULL)
  {
    head= new slist_item<type>(*sl.head);
    tail=head;
    while (tail->next)
    {
      tail= tail->next;
    }
  }
  else
  {
    head= NULL;
    tail= NULL;
  }
}

/***********************************************************************/
/* public member operator=                                             */
/* copies the whole list from the list defined by the argument.        */
/***********************************************************************/

template <class type>
const base_slist<type>& base_slist<type>::operator=(const base_slist<type>& sl)
{
  this->sz= sl.size();
  if (head != NULL)
  {
    delete head;
  }
  if (sl.head != NULL)
  {
    head= new slist_item<type>(*sl.head);
    tail=head;
    while (tail->next)
    {
      tail= tail->next;
    }
  }
  else
  {
    head= NULL;
    tail= NULL;
  }
  return *this;
}

/***********************************************************************/
/* destructor.                                                         */
/***********************************************************************/

template <class type>
base_slist<type>::~base_slist()
{
  slist_item<type> *tmp= head;
  while (head != NULL) {
    tmp= head->next;
    delete head;
    head= tmp;
  }
}

/***********************************************************************/
/* public member function add                                          */
/* prepends an item to the base_slist. It returns a pointer to the     */
/* the prepended slist_item.                                           */
/***********************************************************************/

template <class type>
slist_item<type> *base_slist<type>::add(const type& dat)
{
  if (head != NULL)
  {
    slist_item<type> *it= new slist_item<type>(dat);
    it->next= head;
    head= it;
  }
  else 
  {
    head= new slist_item<type>(dat);
    tail= head;
  }
  this->sz++;
  return head;
}

/***********************************************************************/
/* public mamber function in                                           */
/* takes an argument of type and returns 1 if there is data in the     */
/* base_slist identical to the argument, 0 otherwise.                  */
/***********************************************************************/

template <class type>
int base_slist<type>::in(const type& dat)
{
  slist_item<type> *parent= parent_search(dat);
  if (((parent == NULL) && (this->sz > 0))
      || ((parent) && (parent->next)))
  {
    return 1;
  }
  else 
  {
    return 0;
  }
}

/***********************************************************************/
/* public member function del                                          */
/* takes one argument of type and searches for identical data in the   */
/* base_slist and deletes the item cotaining it if found. It returns   */
/* 1 if anything was deleted, 0 otherwise. Only one item is deleted.   */
/* If there are multiple copies of the specific data in the base_slist */
/* it might remain in the list after del has terminated successfully.  */
/***********************************************************************/

template <class type>
int base_slist<type>::del(const type& dat)
{
  if (head != NULL)
  {
    slist_item<type> *parent= parent_search(dat);
    if (parent)
    {
      if (parent->next)
      {  
	slist_item<type> *it= parent->next;
	parent->next= it->next;
	if (it == tail)
	{
	  tail= parent;
	}
	delete it;
	this->sz--;
	return 1;
      }
    }
    else
    {
      slist_item<type> *it= head;
      head= head->next;
      if (it == tail)
      {
	tail= NULL;
      }
      delete it;
      this->sz--;
      return 1;
    }
  }
  return 0;
}

/***********************************************************************/
/* public member function clear                                        */
/* takes no argument and empties the whole base_slist.                 */
/***********************************************************************/

template <class type>
void base_slist<type>::clear()
{
  this->sz= 0;
  slist_item<type> *tmp= head;
  while (head != NULL) {
    tmp= head->next;
    delete head;
    head= tmp;
  }
  head= NULL;
  tail= NULL;
}

/***********************************************************************/
/* public member function iterator()                                   */
/* returns a freshly created iterator for the base_slist. The iterator */
/* must be deleted by the calling function otherwise it will continue  */
/* to exist as garbage idefinitely.                                    */
/***********************************************************************/

template <class type>
slist_iterator<type> *base_slist<type>::iterator()
{
  return (new slist_iterator<type>(this));
}









