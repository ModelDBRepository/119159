/***********************************************************************/
/* file slset.cc                                                       */
/* contains the implementation of class members of class               */
/* slset.                                                              */
/***********************************************************************/

#include "slset.h"

/***********************************************************************/
/* implementation of class slset                                       */
/***********************************************************************/

/***********************************************************************/
/* public member operator=                                             */
/* copies all data members of the passed slset.                        */
/***********************************************************************/

template <class type>
const slset<type>& slset<type>::operator=(const slset<type>& sl)
{
  this->sz= sl.sz();
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
}


/***********************************************************************/
/* public member function add                                          */
/* takes one argument of type and prepends an item to the slset        */
/* the passed data. If the passed data was already in the set nothing  */
/* changed. add returns a pointer on the prepended slist_item or Null  */
/* if nothing was prepended.                                           */
/***********************************************************************/

template <class type>
slist_item<type> *slset<type>::add(const type& dat)
{
  if (this->head != NULL)
  {
    if (in(dat))
    {
      return NULL;
    }
    else
    {
      slist_item<type> *it= new slist_item<type>(dat);
      it->next= this->head;
      this->head= it;
    }
  }
  else
  {
    this->head= new slist_item<type>(dat);
    this->tail= this->head;
  }
  this->sz++;
  return this->head;
}


/***********************************************************************/
/* public memeber function card                                        */
/* takes no argument and returns the size of the set.                  */
/***********************************************************************/

template <class type>
int slset<type>::card()
{
  return this->sz;
}


/***********************************************************************/
/* public member operator|=                                            */
/* unifies the slset with the passed one. It returns a pointer on the  */
/* owner object itself.                                                */
/***********************************************************************/

template <class type>
slset<type>& slset<type>::operator|=(const slset<type>& s2)
{
  slist_iterator<type> iter(&s2);

  while(iter.step())
  {
    add(iter.c_value());
  }
  return *this;
}


/***********************************************************************/
/* public member operator&=                                            */
/* intersects the slset with the passed slset. It returns a pointer to */
/* the owner object (the slset) itself.                                */
/***********************************************************************/

template <class type>
slset<type>& slset<type>::operator&=(const slset<type>& s2)
{
  slset<type> buf(*this);
  slist_iterator<type> iter(&s2);
  
  clear();
  while(iter.step())
  {
    if (buf.in(iter.c_value()))
    {
      add(iter.c_value());
    }
  }
  return *this;
}


/***********************************************************************/
/* public member operator -=                                           */
/* substracts the passed slset and returns a pointer on the owner      */
/* (the slset) itself.                                                 */
/***********************************************************************/

template <class type>
slset<type>& slset<type>::operator-=(const slset<type>& s2)
{
  slist_iterator<type> iter(&s2); 
  
  while(iter.step())
  {
    del(iter.c_value());
  }
  return *this;
}


/***********************************************************************/
/* public member operator==                                            */
/* tests the equality of th slset with the passed set. slsets are      */
/* considered equal iff they contain exactly the same elements, order  */
/* of appearance in the underlying slist is not considered.            */
/***********************************************************************/


template <class type>
int slset<type>::operator==(const slset<type>& s2)
{
  if (this->sz != s2.sz)
  {
    return 0;
  }
  else
  {
    slist_iterator<type> iter(&s2);
    int equal=1;
    while(iter.step() && equal)
    {
      equal= in(iter.c_value());
    }
    return equal;
  }
}


/***********************************************************************/
/* operator<< puts the slset onto an ostream                           */
/***********************************************************************/


template <class type>
ostream& operator<<(ostream& os, const slset<type>& set)
{
  slist_iterator<type> iter(&set);

  while(iter.step())
  {
    os << iter.c_value() << " ";
  }
  return os;
}

