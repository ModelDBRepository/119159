/***********************************************************************/
/* file bst.cc                                                         */
/* contains the implementation of class members of class               */
/* bst.                                                                */
/***********************************************************************/

#include "bst.h"

/***********************************************************************/
/* implementation of class bst                                         */
/***********************************************************************/

/***********************************************************************/
/* public member operator=                                             */
/* copies the whole bst passed as argument                             */
/***********************************************************************/

template <class type>
const bst<type>& bst<type>::operator=(const bst<type>& tree)
{
  this->sz= tree.size();
  if (this->head != NULL)
  {
    delete this->head;
  }
  if (tree.head)
  {
    this->head= new bst_item<type>(*tree.head);
  }
  else
  {
    this->head= NULL;
  }
  return *this;
}

/***********************************************************************/
/* public member function add                                          */
/* takes one argument of type and inserts an item to the bst,          */
/* containing the passed argument as data. It returns a pointer to the */
/* added bst_item.                                                     */
/***********************************************************************/

template <class type>
bst_item<type> *bst<type>::add(const type& dat)
{
  bst_item<type> *it= NULL;  
  if (this->head != NULL)
  {
    bst_item<type> *parent= parent_search(dat);

    if (parent)
    {
      it= child(parent, dat);
      if (!it)
      {
	it= new bst_item<type>(dat);
	add_item(parent, it);
      }
    }
  }
  else
  {
    it= new bst_item<type>(dat);
    add_item(this->head, it);
  }
  return it;
}


/***********************************************************************/
/* public member function search                                       */
/* takes an argument of type and searches for a bst_item in the bst    */
/* containing it as data. The type type must supply the <>== operators.*/
/* Return value is a pointer to the node containing the argument       */
/* if found, NULL otherwise.                                           */
/***********************************************************************/

template <class type>
bst_item<type> *bst<type>::search(const type& dat)
{
  if (this->head->data == dat)
  {
    return this->head;
  }  
  bst_item<type> *it= NULL;
  bst_item<type> *parent= parent_search(dat);
  if (parent)
  {
    it= child(parent, dat);
  }
  return it;
}

