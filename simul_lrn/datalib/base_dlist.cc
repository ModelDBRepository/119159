/*********************************************************************
 file base_dlist.cc                                                  
 contains the implementation of class members of class               
 base_dlist.                                                         
*********************************************************************/

#include "base_dlist.h"


/*********************************************************************
 implementation of class base_dlist                                  
*********************************************************************/


/*********************************************************************
 private member function find                                        
*********************************************************************/

template <class type>
dlist_item<type> *base_dlist<type>::find(const type &dat)
{
  dlist_item<type> *it;
  it= head;
  while ((it != NULL) && (it->data != dat)) it= it->next;
  return it;
}
  
/*********************************************************************
 constructor                                                         
*********************************************************************/

template <class type>
base_dlist<type>::base_dlist()
{
   head = NULL;
   tail= NULL;
}

/*********************************************************************
 destructor.                                                         
*********************************************************************/

template <class type>
base_dlist<type>::~base_dlist()
{
  clear();
}

/*********************************************************************
 public member function add                                          
 prepends an item to the base_dlist. It returns a pointer to the     
 the prepended dlist_item.                                           
*********************************************************************/

template <class type>
dlist_item<type> *base_dlist<type>::add(const type& dat)
{
  if (head != NULL)
  {
    dlist_item<type> *it= new dlist_item<type>(dat);
    it->next= head;
    head= it;
  }
  else 
  {
    head= new dlist_item<type>(dat);
    tail= head;
  }
  this->sz++;
  return head;
}

/*********************************************************************
 public mamber function in                                           
 takes an argument of type and returns 1 if there is data in the     
 base_dlist identical to the argument, 0 otherwise.                  
*********************************************************************/

template <class type>
int base_dlist<type>::in(const type& dat)
{
  return (find(dat) != NULL);
}

/*********************************************************************
 public member function del_item 
 takes one argument of type and searches for identical data in the   
 base_dlist and deletes the item cotaining it if found. It returns   
 1 if anything was deleted, 0 otherwise. Only one item is deleted.   
 If there are multiple copies of the specific data in the base_dlist 
 it might remain in the list after del has terminated successfully.  
*********************************************************************/

template <class type>
void base_dlist<type>::del_item(dlist_item<type> *it)
{
  assert(it != NULL);
  if (it->prev != NULL) {
    it->prev->next= it->next;
  }
  else {
    head= it->next;
  }
  if (it->next != NULL) {
    it->next->prev= it->prev;
  }
  else {
    tail= it->prev;
  }
  delete it;
  this->sz--;
}


/*********************************************************************
 public member function del                                          
 takes one argument of type and searches for identical data in the   
 base_dlist and deletes the item cotaining it if found. It returns   
 1 if anything was deleted, 0 otherwise. Only one item is deleted.   
 If there are multiple copies of the specific data in the base_dlist 
 it might remain in the list after del has terminated successfully.  
*********************************************************************/

template <class type>
int base_dlist<type>::del(const type& dat)
{
  dlist_item<type> *it= find(dat);
  if ( it != NULL) {
    del_item(it);
    return 1;
  }
  else {
    return 0;
  }
}


/*********************************************************************
 public member function clear                                        
 takes no argument and empties the whole base_dlist.                 
*********************************************************************/

template <class type>
void base_dlist<type>::clear()
{
  this->sz= 0;
  dlist_item<type> *tmp= head;
  while (head != NULL) {
    tmp= head->next;
    delete head;
    head= tmp;
  }
  head= NULL;
  tail= NULL;
}

/*********************************************************************
 public member function iterator()                                   
 returns a freshly created iterator for the base_dlist. The iterator 
 must be deleted by the calling function otherwise it will continue  
 to exist as garbage idefinitely.                                    
*********************************************************************/

template <class type>
dlist_iterator<type> *base_dlist<type>::iterator()
{
  return (new dlist_iterator<type>(this));
}









