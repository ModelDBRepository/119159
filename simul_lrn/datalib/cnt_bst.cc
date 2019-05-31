/***********************************************************************/
/* file cnt_bst.cc                                                     */
/* contains the implementation of class members of class               */
/* cnt_bst.                                                            */
/***********************************************************************/

#include "cnt_bst.h"

/***********************************************************************/
/* implementation of class cnt_bst                                     */
/***********************************************************************/

/***********************************************************************/
/* public member operator=                                             */
/* copies all data members' contents of the passed cnt_bst             */
/***********************************************************************/

template <class type>
const cnt_bst<type>& cnt_bst<type>::operator=(const cnt_bst<type>& tree)
{
  this->sz= tree.sz;
  if (this->head != NULL)
  {
    delete this->head;
  }
  if (tree.head)
  {
    this->head= new cnt_bst_item<type>(*tree.head);
  }
  else
  {
    this->head= NULL;
  }
}

/***********************************************************************/
/* public member function add                                          */
/* takes one argument of type and inserts an item to the bst           */
/* if it doesn't already exist. Otherwise the count of the existing    */
/* item is increased. It returns a pointer to the added bst_item.      */
/***********************************************************************/

template <class type>
cnt_bst_item<type> *cnt_bst<type>::add(const type& dat)
{
  return cnt_add(dat, 1);
}

/***********************************************************************/
/* public member function cnt_add                                      */
/* works like add (in doubleity it's vice versa ...) but takes an        */
/* argument of int which denotes the count for the added item (or the  */
/* increase of count if it already existed in the tree).               */
/* As is easy to see add(dat) is simply cnt_add(dat,1).                */
/***********************************************************************/

template <class type>
cnt_bst_item<type> *cnt_bst<type>::cnt_add(const type& dat, int cnt)
{
  cnt_bst_item<type>* cnt_it= NULL;
  if (this->head != NULL)
  {
    bst_item<type> *parent= parent_search(dat);
    
    if (parent)
    {
      bst_item<type> *it= child(parent, dat);
      if (it)
      {
	cnt_it= (cnt_bst_item<type> *) it;
	cnt_it->count+= cnt;
      }
      else
      {
	cnt_it= new cnt_bst_item<type>(dat, cnt);
	add_item(parent, cnt_it);
      }
    }
    else
    {
      cnt_it= (cnt_bst_item<type> *) this->head;
      cnt_it->count+= cnt;
    }   
  }
  else
  {
    cnt_it= new cnt_bst_item<type>(dat, cnt);
    add_item(this->head, cnt_it);
  }
  return cnt_it;
}

/***********************************************************************/
/* public member function search                                       */
/* takes an argument of type and returns a pointer to the node where   */
/* the first argument was found containing the passed argument as      */
/* data, NULL if it wasn't found. type must supply == < > operators.   */
/***********************************************************************/

template <class type>
cnt_bst_item<type> *cnt_bst<type>::search(const type& dat)
{
  if (this->head->data == dat)
  {
    return (cnt_bst_item<type> *) this->head;
  }
  cnt_bst_item<type> *it= NULL;
  bst_item<type> *parent= parent_search(dat);
  if (parent)
  {
    it= (cnt_bst_item<type> *) child(parent, dat);
  }
  return it;
}

/***********************************************************************/
/* protected member function del_item                                  */
/* workhorse for del, deletes an item                                  */
/* overloads the according function of base_bst                        */
/***********************************************************************/

template <class type>
void cnt_bst<type>::del_item(bst_item<type> *parent,
			      bst_item<type> *it) 
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
    ((cnt_bst_item<type> *) it)->count=
	((cnt_bst_item<type> *) next_it)->count;
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
      ((cnt_bst_item<type> *) it)->count=
	  ((cnt_bst_item<type> *) next_it)->count;
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
      if (it == this->head)
      {
	this->head= NULL;
      }
      delete it;
      this->sz--;
    }
  }
}

/***********************************************************************/
/* public member function del                                          */
/* takes a argument of type and deletes the node containing identical  */
/* data to the passed argument. It returns 1 if anything was deleted,  */
/* 0 otherwise. type must supply == < > operators.                     */
/***********************************************************************/

template <class type>
int cnt_bst<type>::del(const type& dat, int cnt)
{
  if (this->head != NULL)
  {
    bst_item<type> *parent= parent_search(dat);
    if (parent)
    {
      bst_item<type> *it= child(parent,dat);   
      if (it)
      {
	int new_cnt= ((cnt_bst_item<type> *) it)->count - cnt;
	if (new_cnt < 1)
	{
	  del_item(parent, it);
	}
	else
	{
	  ((cnt_bst_item<type> *) it)->count= new_cnt;
	}
	return 1;
      }
    }
    else
    {
      bst_item<type>* it= this->head;
      int new_cnt= ((cnt_bst_item<type> *) it)->count - cnt;
      if (new_cnt < 1)
      {
	del_item(parent, it);
      }
      else
      {
	((cnt_bst_item<type> *) it)->count= new_cnt;
      }
      return 1;
    }
  }
  return 0;
}

/***********************************************************************/
/* public member function count                                        */
/* takes one argument of type and returns the count it has in the      */
/* cnt_bst, 0 if it isn't in the cnt_bst.                              */
/***********************************************************************/

template <class type>
int cnt_bst<type>::count(const type &dat)
{
  cnt_bst_item<type> *it= search(dat);
  if (it)
  {
    return it->count;
  }
  else
  {
    return 0;
  }
}

/***********************************************************************/
/* public member function iterator()                                   */
/* returns a pointer to a freshly created iterator for the cnt_bst.    */
/* The Iterator must be deleted by the calling function or it will     */
/* continue to exist as garbage indefinitely.                          */
/***********************************************************************/

template <class type>
cnt_bst_iterator<type> *cnt_bst<type>::iterator()
{
  return (new cnt_bst_iterator<type>(this));
}

