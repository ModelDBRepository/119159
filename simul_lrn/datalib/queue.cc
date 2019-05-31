/***********************************************************************/
/* file queue.cc                                                       */
/* contains the implementation of class members of class               */
/* queue.                                                              */
/***********************************************************************/

#include "queue.h"

/***********************************************************************/
/* implementation of class queue                                       */
/***********************************************************************/

/***********************************************************************/
/* operator=                                                           */
/* copies all data members of the passed queue                         */
/***********************************************************************/

template <class type>
const queue<type>& queue<type>::operator=(const queue<type>& qu)
{
  slist_item<type> *tmp;
  if (this->head != NULL)
  {
    this->clear();
  }
  this->sz= qu.sz;
  if (qu.head)
  {
    tmp= qu.head;
    do {
      push(tmp->data);
      tmp= tmp->next;
    } while (tmp);
  }
  else
  {
    this->head= NULL;
    this->tail= NULL;
  }
  return *this;
}


/***********************************************************************/
/* public member function push                                         */
/* takes an argument of type and adds a new item containing the data   */
/* passed to the top of the queue                                      */
/***********************************************************************/

template <class type>
void queue<type>::push(const type& dat)
{
  if (this->head != NULL)           // same as if (this->tail != NULL) 
  {
    this->tail->next= new slist_item<type>(dat);
    this->tail= this->tail->next;
  }
  else
  {
    this->head= new slist_item<type>(dat);
    this->tail= this->head;
  }
  this->sz++;
}


/***********************************************************************/
/* public mamber function pop                                          */
/* takes no argument and returns the data of the last item in the      */
/* queue. The item is deleted                                          */
/***********************************************************************/

template <class type>
type queue<type>::pop()
{
  assert(this->head != NULL);
  type dat= this->head->data;
  slist_item<type> *it= this->head->next;
  this->head->next= NULL;
  delete this->head;
  this->head=it;
  this->sz--;
  return dat;
}


/***********************************************************************/
/* public member function peek                                         */
/* returns the data contained in the first item in the queue.          */
/***********************************************************************/

template <class type>
type queue<type>::peek()
{
  assert(this->head != NULL);
  return this->head->data;
}

/***********************************************************************/
/* public member function peek                                         */
/* returns the data contained in the first item in the queue.          */
/***********************************************************************/

template <class type>
type queue<type>::last()
{
  assert(this->head != NULL);
  return this->tail->data;
}

