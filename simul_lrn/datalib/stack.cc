/***********************************************************************/
/* file stack.cc                                                       */
/* contains the implementation of class members of class               */
/* stack.                                                              */
/***********************************************************************/

#include "stack.h"

/***********************************************************************/
/* implementation of class stack                                       */
/***********************************************************************/

/***********************************************************************/
/* public member operator=                                             */
/* copies all data members of the passed stack.                        */
/***********************************************************************/

template <class type>
const stack<type>& stack<type>::operator=(const stack<type>& st)
{
  this->sz= st.sz;
  if (this->head != NULL)
  {
    delete this->head;
  }
  if (st.head)
  {
    this->head= new slist_item<type>(*st.head);
    this->tail=this->head;
    while (this->tail->next)
    {
      this->tail= this->tail->next;
    }
  }
  else
  {
    this->head= NULL;
    this->tail= NULL;
  }
  return *this;
}


/***********************************************************************/
/* public member funtion push                                          */
/* takes one argument of type and prepends an item with according      */
/* data using the more basic member function add.                      */ 
/***********************************************************************/

template <class type>
void stack<type>::push(const type& dat)
{
  add(dat);
}


/***********************************************************************/
/* public member function pop                                          */
/* takes no argument and returns the content of the topmost item in    */
/* the stack. The topmost item is deleted.                             */
/***********************************************************************/

template <class type>
type stack<type>::pop()
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
/* returns the content of the topmost item in the stack without        */
/* deleting it.                                                        */
/***********************************************************************/

template <class type>
type stack<type>::peek()
{
  assert(this->head != NULL);
  return this->head->data;
}












