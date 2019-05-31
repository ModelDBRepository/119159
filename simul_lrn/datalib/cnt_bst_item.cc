/***********************************************************************/
/* file cnt_bst_item.cc                                                */
/* contains the implementation of class members of class               */
/* cnt_bst_item.                                                       */
/***********************************************************************/

#include "cnt_bst_item.h"

/***********************************************************************/
/* implementation of class cnt_bst_item                                */
/***********************************************************************/

/***********************************************************************/
/* constructor                                                         */
/***********************************************************************/

template <class type>
cnt_bst_item<type>::cnt_bst_item(const type& dat, int cnt):
  bst_item<type>(dat)
{
  count= cnt;
}

/***********************************************************************/
/* copy constructor                                                    */
/***********************************************************************/

template <class type>
cnt_bst_item<type>::cnt_bst_item(const cnt_bst_item<type>& it):
  bst_item<type>(it)
{
  count= it.count;
}

/***********************************************************************/
/* public member operator=                                             */
/* copies the data members of the passed cnt_bst_item.                 */
/***********************************************************************/

template <class type>
const cnt_bst_item<type>& cnt_bst_item<type>::operator=
(const cnt_bst_item<type>& it)
{
  this->data= it.data;
  count= it.count;
  if (this->left)
  {
    delete this->left;
  }
  if (it.left)
  {
    this->left= new bst_item<type>(*it.left);
  }
  else
  {
    this->left= NULL;
  }
  if (this->right)
  {
    delete this->right;
  }
  if (it.right)
  {
    this->right= new bst_item<type>(*it.right);
  }
  else
  {
    this->right= NULL;
  }
  return it;
}
 


