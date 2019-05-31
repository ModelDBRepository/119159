/*********************************************************************
 file base_label_dlist.cc                                                  
 contains the implementation of class members of class               
 base_label_dlist.                                                         
*********************************************************************/

#include "base_label_dlist.h"


/*********************************************************************
 implementation of class base_label_dlist                                  
*********************************************************************/


/*********************************************************************
 private member function find                                        
*********************************************************************/

template <class type>
label_dlist_item<type> *base_label_dlist<type>::find(const type &dat)
{
  label_dlist_item<type> *it;
  it= head;
  while ((it != NULL) && (it->data != dat)) it= it->next;
  return it;
}
  
/*********************************************************************
 private member function find                                        
*********************************************************************/

template <class type>
label_dlist_item<type> *base_label_dlist<type>::label_find(int lab)
{
  label_dlist_item<type> *it;
  it= head;
  while ((it != NULL) && (it->label < lab)) it= it->next;
  if ((it) && (it->label == lab)) return it;
  else return NULL;
}

/*********************************************************************
 private member function find                                        
*********************************************************************/

template <class type>
type base_label_dlist<type>::label_find_data(int lab)
{
  label_dlist_item<type> *it;
  it= head;
  while ((it != NULL) && (it->label < lab)) it= it->next;
  if ((it) && (it->label == lab)) return it->data;
  else return NULL;
}
  
/*********************************************************************
 constructor                                                         
*********************************************************************/

template <class type>
base_label_dlist<type>::base_label_dlist()
{
   head = NULL;
   tail= NULL;
}

/*********************************************************************
 destructor.                                                         
*********************************************************************/

template <class type>
base_label_dlist<type>::~base_label_dlist()
{
  dlist_item<type> *tmp= head;
  while (head != NULL) {
    tmp= head->next;
    delete head;
    head= tmp;
  }
}


/*********************************************************************
 public member function insert                                       
 takes a pointer to a label_dlist node and an argument type. It inserts a  
 new item with the correct data behind the node received as first    
 argument.                                                            
*********************************************************************/

template <class type>
label_dlist_item<type> *base_label_dlist<type>::insert(label_dlist_item<type> *parent, int label, const type& dat)
{
  assert(parent != NULL);
  label_dlist_item<type> *it= new label_dlist_item<type>(label, dat);
  it->next= parent->next;
  it->prev= parent;
  parent->next= it;
  if (it->next != NULL) it->next->prev= it;
  else tail= it;
  this->sz++;

  return it;
}

/*********************************************************************
 public member function prepend
 takes a label and argument of type. It inserts a  
 new item with the correct data at the beginning of the list                  
*********************************************************************/

template <class type>
label_dlist_item<type> *base_label_dlist<type>::prepend(int label, const type & dat)
{
  head->prev= new label_dlist_item<type>(label, dat);
  head->prev->next= head;
  head= head->prev;
  head->prev= NULL;
  return head;
}

/*********************************************************************
 public member function add                                          
 prepends an item to the base_label_dlist. It returns a pointer to the     
 the prepended label_dlist_item.                                           
*********************************************************************/

template <class type>
label_dlist_item<type> *base_label_dlist<type>::add(int label, const type& dat)
{
  if (head != NULL) {
    label_dlist_item<type> *prnt= head;
    while ((prnt->next != NULL) && (prnt->next->label < label))
      prnt= prnt->next;
    if (prnt->label < label) {
      return insert(prnt, label, dat);
    }
    else { // we only have larger stuff in the list so far ...
      return prepend(label, dat);
    }
  }
  else {
    head= new label_dlist_item<type>(label, dat);
    tail= head;
    this->sz++;
    return head;
  }
}
      
/*********************************************************************
 public mamber function in                                           
 takes an argument of type and returns 1 if there is data in the     
 base_label_dlist identical to the argument, 0 otherwise.                  
*********************************************************************/

template <class type>
int base_label_dlist<type>::in(const type& dat)
{
  return (find(dat) != NULL);
}

/*********************************************************************
 public mamber function label in                                           
 takes an argument of type and returns 1 if there is data in the     
 base_label_dlist identical to the argument, 0 otherwise.                  
*********************************************************************/

template <class type>
int base_label_dlist<type>::label_in(int label)
{
  return (label_find(label) != NULL);
}

/*********************************************************************
 public member function del_item 
 takes one argument of type and searches for identical data in the   
 base_label_dlist and deletes the item cotaining it if found. It returns   
 1 if anything was deleted, 0 otherwise. Only one item is deleted.   
 If there are multiple copies of the specific data in the base_label_dlist 
 it might remain in the list after del has terminated successfully.  
*********************************************************************/

template <class type>
void base_label_dlist<type>::del_item(label_dlist_item<type> *it)
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
 base_label_dlist and deletes the item cotaining it if found. It returns   
 1 if anything was deleted, 0 otherwise. Only one item is deleted.   
 If there are multiple copies of the specific data in the base_label_dlist 
 it might remain in the list after del has terminated successfully.  
*********************************************************************/

template <class type>
int base_label_dlist<type>::del(const type& dat)
{
  label_dlist_item<type> *it= find(dat);
  if ( it != NULL) {
    del_item(it);
    return 1;
  }
  else {
    return 0;
  }
}


/*********************************************************************
 public member function del                                          
 takes one argument of type and searches for identical data in the   
 base_label_dlist and deletes the item cotaining it if found. It returns   
 1 if anything was deleted, 0 otherwise. Only one item is deleted.   
 If there are multiple copies of the specific data in the base_label_dlist 
 it might remain in the list after del has terminated successfully.  
*********************************************************************/

template <class type>
int base_label_dlist<type>::label_del(int label)
{
  label_dlist_item<type> *it= label_find(label);
  if (it != NULL) {
    del_item(it);
    return 1;
  }
  else {
    return 0;
  }
}

/*********************************************************************
 public member function clear                                        
 takes no argument and empties the whole base_label_dlist.                 
*********************************************************************/

template <class type>
void base_label_dlist<type>::clear()
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
 returns a freshly created iterator for the base_label_dlist. The iterator 
 must be deleted by the calling function otherwise it will continue  
 to exist as garbage idefinitely.                                    
*********************************************************************/

template <class type>
label_dlist_iterator<type> *base_label_dlist<type>::iterator()
{
  return (new label_dlist_iterator<type>(this));
}









