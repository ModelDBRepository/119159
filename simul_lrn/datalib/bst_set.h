/***********************************************************************/
/* file bst_set.h                                                      */
/* contains the class definition of class bst_set.                     */
/***********************************************************************/

#ifndef BST_SET_H
#define BST_SET_H

#include <assert.h>
#include "bst.h"

/***********************************************************************/
/* class definition of class bst_set. bst_set is a set class based on  */
/* a bst.                                                              */
/***********************************************************************/

template <class type>
class bst_set: public bst<type>
{
   friend class bst_iterator<type>;              // bst_iterator may access
                                                 // private data fields
public:
   bst_set() { }                                 // constructor
   bst_set(const bst_set<type>& tree):           // copy constructor
     bst<type>(tree) { }
   const bst_set<type>& operator=                // copy operator
     (const bst_set<type>&);
   ~bst_set() { }                                // destructor
   int card();                                   // returns the size
   type min();                                   // return minimum value
   type max();                                   // return maximum value
   type mean();                                  // return the mean
   const bst_set<type>& operator|=               // unification
     (const bst_set<type>&);
   const bst_set<type>& operator&=               // intersection
     (const bst_set<type>&);
   const bst_set<type>& operator-=               // substraction
     (const bst_set<type>&);
   int operator==(const bst_set<type>&);         // tests equality
};

template <class type>
istream &operator>>(istream &, bst_set<type> &);

template <class type>
ostream &operator<<(ostream &, bst_set<type> &);

#include "bst_set.cc"                            // because of templates

#endif

