/***********************************************************************/
/* file bst_item.cc                                                    */
/* contains the implementation of class members of class               */
/* bst_item.                                                           */
/***********************************************************************/

#include "bst_item.h"

/***********************************************************************/
/* implementation of class bst_item                                    */
/***********************************************************************/

/***********************************************************************/
/* constructor                                                         */
/***********************************************************************/

template <class type>
bst_item<type>::bst_item(const type& dat):container_item<type>(dat)
{
   left= NULL;
   right= NULL;
}

/***********************************************************************/
/* copy constructor                                                    */
/* creates also copies of the children recursively                     */
/***********************************************************************/

template <class type>
bst_item<type>::bst_item(const bst_item<type>& it):
  container_item<type>(it.data)
{
  if (it.left)
  {
    left= new bst_item<type>(*it.left);
  }
  else
  {
    left= NULL;
  }
  if (it.right)
  {
    right= new bst_item<type>(*it.right);
  } 
  else
  {
    right= NULL;
  }
}

/***********************************************************************/
/* public member operator=                                             */
/* copies the arguments data into the bst_item and generates children  */
/* identical to those of the argument. The process repeats recursively.*/
/* Prior children (subtrees) of bst_item are deleted recursively.      */
/***********************************************************************/

template <class type>
const bst_item<type>& bst_item<type>::operator=
(const bst_item<type>& it)
{
  this->data= it.data;
  if (left)
  {
    delete left;
  }
  if (it.left)
  {
    left= new bst_item<type>(*it.left);
  }
  else
  {
    left= NULL;
  }
  if (right)
  {
    delete right;
  }
  if (it.right)
  {
    right= new bst_item<type>(*it.right);
  }
  else
  {
    right= NULL;
  }
  return *this;
}


/***********************************************************************/
/* destructor                                                          */
/* destroys all children recursively.                                  */
/***********************************************************************/

template <class type>
bst_item<type>::~bst_item()
{
  if (left)
  {
    delete left;
  }
  if (right)
  {
    delete right;
  }
}




