/***********************************************************************/
/* file cnt_bst_set.cc                                                 */
/* contains the implementation of class members of class               */
/* cnt_bst_set.                                                        */
/***********************************************************************/

#include "cnt_bst_set.h"

/***********************************************************************/
/* implementation of class cnt_bst_set                                 */
/***********************************************************************/

/***********************************************************************/
/* public member operator=                                             */
/* copies all contents of the passed set                               */
/***********************************************************************/

template <class type>
const cnt_bst_set<type>& cnt_bst_set<type>::operator=
(const cnt_bst_set<type>& tree)
{
  this->sz= tree.size();
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
/* public member function card                                         */
/* takes no argument and returns the cardinality (size) of the set.    */
/***********************************************************************/

template <class type>
int cnt_bst_set<type>::card()
{
  return this->sz;
}

/***********************************************************************/
/* public member function max                                          */
/* takes no argument and returns the maximum element of the set.       */
/***********************************************************************/

template <class type>
type cnt_bst_set<type>::max()
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
/* public member function minimum                                      */
/* takes no argument and returns the minimum element of the set.       */
/***********************************************************************/

template <class type>
type cnt_bst_set<type>::min()
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
/* takes no argument and returns the (weighted) mean of the set.       */
/***********************************************************************/

template <class type>
type cnt_bst_set<type>::mean()
{
  assert(this->head != NULL);
  cnt_bst_iterator<type> iter(this);
  type avg= 0;
  int no= 0;
  while(iter.step())
  {
    avg+= iter.c_value() * iter.c_count();
    no+= iter.c_count();
  }
  avg/= no;
  return avg;
}

/***********************************************************************/
/* public member operator|=                                            */
/* unifies the set with the passed set. If elements exist in both sets */
/* simultaneously, the counts are added.                               */
/***********************************************************************/

template <class type>
cnt_bst_set<type>& cnt_bst_set<type>::operator|=(const cnt_bst_set<type>& s2)
{
  cnt_bst_iterator<type> iter(&s2);

  while(iter.step())
  {
    cnt_add(iter.c_value(), iter.c_count());
  }
  return *this;
}

/***********************************************************************/
/* public member operator&=                                            */
/* intersects the set with the passed set. It is not very well defined */
/* yet, because elements ar all only added once neglecting their count.*/
/***********************************************************************/

template <class type>
cnt_bst_set<type>& cnt_bst_set<type>::operator&=(const cnt_bst_set<type>& s2)
{
  cnt_bst_set<type> buf(*this);
  cnt_bst_iterator<type> iter(&s2);
  
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
/* takes the difference set \ passed set. This means that the count of */
/* all elements in the set are diminished by their count in the passed */
/* set and deleted completely if this exceeds their current count.     */
/***********************************************************************/

template <class type>
cnt_bst_set<type>& cnt_bst_set<type>::operator-=(const cnt_bst_set<type>& s2)
{
  cnt_bst_iterator<type> iter(&s2);
  
  while(iter.step())
  {
    del(iter.c_value(), iter.c_count());
  }
  return *this;
}

/***********************************************************************/
/* public member operator==                                            */
/* tests equality to the passed set. Sets are considered equal iff     */
/* they have the same elements. This is no good definition yet, should */
/* considered equal iff elements and counts coincide...be              */
/***********************************************************************/


template <class type>
int cnt_bst_set<type>::operator==(const cnt_bst_set<type>& s2)
{
  if (this->sz != s2.sz)
  {
    return 0;
  }
  else
  {
    cnt_bst_iterator<type> iter(&s2);
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
istream& operator>>(istream& is, cnt_bst_set<type>& set)
{
  int nn;
  type tmp;

  while (is.good())
  {
    is >> tmp;
    if (is.good())
    {
      is >> nn;
      set.cnt_add(tmp, nn);
    }
  }

  return is;
}


/***********************************************************************/
/* operator<< puts the set onto an ostream                             */
/***********************************************************************/


template <class type>
ostream& operator<<(ostream& os, const cnt_bst_set<type>& set)
{
  cnt_bst_iterator<type> iter(&set);

  while(iter.step())
  {
    os << iter.c_value() << " " << iter.c_count() << endl;
  }
  os << endl;

  return os;
}

