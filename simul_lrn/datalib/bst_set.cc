/***********************************************************************/
/* file bst_set.cc                                                     */
/* contains the implementation of class members of class               */
/* bst_set.                                                            */
/***********************************************************************/

#include "bst_set.h"

/***********************************************************************/
/* implementation of class bst_set                                     */
/***********************************************************************/

/***********************************************************************/
/* public member operator=                                             */
/* creates a new set with identical content as the passed argument.    */
/***********************************************************************/

template <class type>
const bst_set<type>& bst_set<type>::operator=(const bst_set<type>& tree)
{
  this->sz= tree.sz;
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
/* public member function card                                         */
/* takes no argument and returns the cardinality (size) of the set.    */
/***********************************************************************/

template <class type>
int bst_set<type>::card()
{
  return this->sz;
}

/***********************************************************************/
/* public member function max                                          */
/* takes no argument and returns the maximum value contained in the    */
/* set, i.e. the data of the rightmost item in the underlying tree.    */
/***********************************************************************/

template <class type>
type bst_set<type>::max()
{
  assert(this->head != NULL);
  bst_item<type> *it= this->head;
  while (it->right)
  {
    it= it->right;
  }
  return it->data;
}

/***********************************************************************/
/* public member function min                                          */
/* takes no argument and returns the minimum value contained in the    */
/* set, i.e. the data of the leftmost item in the underlying tree.     */
/***********************************************************************/

template <class type>
type bst_set<type>::min()
{
  assert(this->head != NULL);
  bst_item<type> *it= this->head;
  while (it->left)
  {
    it= it->left;
  }
  return it->data;
}


/***********************************************************************/
/* public member function mean                                         */
/* takes no argument and returns the mean of the set.                  */
/***********************************************************************/

template <class type>
type bst_set<type>::mean()
{
  assert(this->head != NULL);
  bst_iterator<type> iter(this);
  type avg= 0;
  while(iter.step())
  {
    avg+= iter.c_value();
  }
  avg/= this->sz;
  return avg;
}



/***********************************************************************/
/* public member operator|=                                            */
/* unites the set with the argument. The results is stored in the set  */
/* itself.                                                             */
/***********************************************************************/

template <class type>
const bst_set<type>& bst_set<type>::operator|=(const bst_set<type>& s2)
{
  bst_iterator<type> iter(&s2);

  while(iter.step())
  {
    add(iter.c_value());
  }
  return *this;
}


/***********************************************************************/
/* public member operator&=                                            */
/* intersects with the argument and stores the result in the set       */
/* itself.                                                             */
/***********************************************************************/

template <class type>
const bst_set<type>& bst_set<type>::operator&=(const bst_set<type>& s2)
{
  bst_set<type> buf(*this);
  bst_iterator<type> iter(&s2);
  
  this->clear();
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
/* substracts the argument and stores the result in the set itself.    */
/***********************************************************************/

template <class type>
const bst_set<type>& bst_set<type>::operator-=(const bst_set<type>& s2)
{
  bst_iterator<type> iter(&s2); 
  
  while(iter.step())
  {
    this->del(iter.c_value());
  }
  return *this;
}

/***********************************************************************/
/* public member operator==                                            */
/* tests equality to the argument set. Set are considered equal iff    */
/* they have exactly the same elements (not necessarily in the same    */
/* order in the underlying tree).                                      */
/***********************************************************************/

template <class type>
int bst_set<type>::operator==(const bst_set<type>& s2)
{
  if (this->sz != s2.sz)
  {
    return 0;
  }
  else
  {
    bst_iterator<type> iter(&s2);
    int equal=1;
    while(iter.step() && equal)
    {
      equal= in(iter.c_value());
    }
    return equal;
  }
}

/***********************************************************************/
/* operator>>                                                          */
/* gets the set from the passed istream.                               */
/***********************************************************************/

template <class type>
istream& operator>>(istream& is, bst_set<type>& set)
{
  int nn, i;
  type tmp;

  is >> nn;
  for (i= 0; i < nn; i++)
  {
    is >> tmp;
    set.add(tmp);
  }

  return is;
}



/***********************************************************************/
/* operator<<                                                          */
/* outputs the set onto the passed ostream.                            */
/***********************************************************************/

template <class type>
ostream& operator<<(ostream& os, bst_set<type>& set)
{
  bst_iterator<type> iter(&set);

  os << set.size() << endl;
  while(iter.step())
  {
    os << iter.c_value() << " ";
  }
  os << endl;
  return os;
}

