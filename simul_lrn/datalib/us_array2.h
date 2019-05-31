/*********************************************************************
 file us-array.h                                                        
 contains the class definition of class us_array (ultra sparse array)
 which is build with label_dlists
*********************************************************************/

#ifndef US_ARRAY_H
#define US_ARRAY_H

#include <assert.h>

template <class type1, class type2>
class tupel
{
 public:
  type1 x;
  type2 y;
};

template <class type1, class type2>
ostream &operator<<(ostream &os, tupel<type1,type2> tp)
{
  os << tp.x << " ; ";
  os << tp.y << " . ";
}

template <class type1, class type2>
int operator==(tupel<type1,type2> tp1, tupel<type1,type2> tp2)
{
  return ((tp1.x == tp2.x) && (tp1.y == tp2.y));
}
	      
template <class type>
class us_array;                   // forward declaration

#include "label_dlist.h"
#include "slset.h"

/*********************************************************************
 class definition of class us_array. us_array is an iterated doubly
 linked list. (us_array = ultra sparse array) 
*********************************************************************/

template <class type>
class us_array
{
  //   friend class us_array_iterator<type>;  // us_array_iterator may access
                                              // private data fields
 protected:
  int dim;
  int sz;
  label_dlist<void *> *thelist;
  
public:
   us_array(int);                                    // constructor
   virtual ~us_array();                              // destructor
   virtual void add(int *, const type &); // adds
   virtual int in(int *);
   virtual int del(int *);          // delete
   virtual void neigh_add(int *, int *, int, label_dlist<type> *, slset<tupel<int *, type> > *); 
   virtual void neigh(int *, int *, int, label_dlist<void *> *, int, slset<tupel<int *, type> > *); // returns set of neighbors
   virtual void neighbors(int *, int, slset<tupel<int*, type> > *);
};

#include "us_array2.cc"                            // because of templates

#endif
