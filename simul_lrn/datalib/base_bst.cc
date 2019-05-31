/***********************************************************************/
/* file base_bst.cc                                                    */
/* contains the implementation of class members of class               */
/* base_bst.                                                           */
/***********************************************************************/

#include "base_bst.h"

/***********************************************************************/
/* implementation of class base_bst                                    */
/***********************************************************************/

/***********************************************************************/
/* constructor                                                         */
/***********************************************************************/

template <class type>
base_bst<type>::base_bst(const base_bst<type>& tree):container<type>(tree)
{
  if (tree.head)
  {
    head= new bst_item<type>(*tree.head);
  }
  else
  {
    head= NULL;
  }
}

/***********************************************************************/
/* public operator=                                                    */
/* copies its arguments' data members on the data members of           */
/* its owner.                                                          */
/***********************************************************************/ 

template <class type>
const base_bst<type>& base_bst<type>::operator=(const base_bst<type>& tree)
{
  this->sz= tree.size();
  if (head != NULL)
  {
    delete head;
  }
  if (tree.head)
  {
    head= new bst_item<type>(*tree.head);
  }
  else
  {
    head= NULL;
  }
  return *this;
}

/***********************************************************************/
/* protected member function del_item                                  */
/* workhorse for del, deletes an item                                  */
/***********************************************************************/

template <class type>
void base_bst<type>::del_item(bst_item<type> *parent, bst_item<type> *it)
{
  assert(it != NULL);
  if (it->left)
  {
    bst_item<type> *old_it= it;
    bst_item<type> *next_it= it->left;
    while(next_it->right)
    {
      old_it= next_it;
      next_it= next_it->right;
    }
    it->data= next_it->data;
    del_item(old_it, next_it);
  }
  else
  {
    if (it->right)
    {
      bst_item<type> *old_it= it;
      bst_item<type> *next_it= it->right;
      while(next_it->left)
      {
	old_it= next_it;
	next_it= next_it->left;
      }
      it->data= next_it->data;
      del_item(old_it, next_it);
    }
    else
    {
      if (parent)
      {
	if (it == parent->left)
	{
	  parent->left= NULL;
	}
	else
	{
	  parent->right= NULL;
	}
      }
      if (it == head)
      {
	head= NULL;
      }
      delete it;
      this->sz--;
    }
  }
}

/***********************************************************************/
/* protected member function internal search takes an argument of type */
/* which specifies what data to search for and an argument of int      */
/* in which 0 for 'not found' and 1 for 'found' is return. The return  */
/* value is either a pointer on the node where the data was found or if*/
/* it was not found a pointer on the node where the search ended.      */
/***********************************************************************/

template <class type>
bst_item<type> *base_bst<type>::parent_search(const type& dat)
{
  bst_item<type> *it= head;
  bst_item<type> *parent= NULL;

  while (it)
  {
    if (it->data == dat)
    {
      return parent;
    }
    else
    {
      parent= it;
      if (dat < it->data)
      {
	it= it->left;
      }
      else
      {
	it= it->right;
      }
    }
  }
  return parent;
}

/***********************************************************************/
/* protected member function child                                     */
/* takes a pointer to bst_item and an argument of type. It determines  */
/* whether the bst_item denoted by argument 1 has a child with data    */
/* identical to argument 2. If so it returns a pointer to the child    */
/* otherwise NULL.                                                     */
/***********************************************************************/

template <class type>
bst_item<type> *base_bst<type>::child(bst_item<type> *parent, const type& dat)
{
  assert(parent != NULL);
  if ((parent->left) && (dat == parent->left->data))
  {
    return parent->left;
  }
  else
  {
    if ((parent->right) && (dat == parent->right->data))
    {
      return parent->right;
    }
    else
    {
      return NULL;
    }
  }
}


/***********************************************************************/
/* protected member function add_item                                  */
/* takes two pointer to bst_item and inserts the item pointed to by    */
/* argument 2 as a child of the item pointed to by argument 1.         */
/* There is no check if there already is a child.                      */
/***********************************************************************/

template <class type>
void base_bst<type>::add_item(bst_item<type> *parent, bst_item<type> *it)
{
  if (it)
  {
    this->sz++;
    if (parent)
    {
      if (it->data < parent->data)
      {
	parent->left= it;
      }
      else
      {
	parent->right= it;
      }      
    }
    else
    {
      head= it;
    }
  }
}


/***********************************************************************/
/* constructor                                                         */
/***********************************************************************/

template <class type>
base_bst<type>::base_bst()
{
  head= NULL;
}


/***********************************************************************/
/* destructor                                                          */
/***********************************************************************/

template <class type>
base_bst<type>::~base_bst()
{
  if (head != NULL) 
  {
    delete head;
  }
}

/***********************************************************************/
/* public member function in                                           */
/* returns 1 if the data passed as single argument is found int the    */
/* base_bst, 0 otherwise.                                              */
/***********************************************************************/

template <class type>
int base_bst<type>::in(const type& dat)
{
  bst_item<type>* parent= parent_search(dat);

  if (((parent == NULL) && (this->sz > 0))
      || ((parent) && (child(parent,dat))))
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
/* takes an argument of type and deletes the node in the base_bst      */
/* containing identical data as the argument. It returns 1 if anything */
/* was deleted, 0 otherwise.                                           */
/***********************************************************************/

template <class type>
int base_bst<type>::del(const type& dat)
{
  if (head != NULL)
  {
    bst_item<type> *parent= parent_search(dat);
    if (parent)
    {
      bst_item<type> *it= child(parent,dat);   
      if (it)
      {
	del_item(parent, it);
	return 1;
      }
    }
    else
    {
      del_item(parent, head);
      return 1;
    }
  }
  return 0;
}

/***********************************************************************/
/* public member function clear                                        */
/* takes no arguments and  empties the whole base_bst                  */
/***********************************************************************/

template <class type>
void base_bst<type>::clear()                       
{
  this->sz= 0;
  if (head != NULL)
  {
    delete head;
  }
  head= NULL;
}

/***********************************************************************/
/* public member function iterator()                                   */
/* takes no argument and returns a freshly created iterator for the    */
/* base_bst. The iterator must be deleted by the calling function.     */
/* Otherwise it will continue as garbage indefinitely.                 */
/***********************************************************************/

template <class type>
bst_iterator<type> *base_bst<type>::iterator()
{
  bst_iterator<type> *tmp= new bst_iterator<type>(this);
  return tmp;
}
