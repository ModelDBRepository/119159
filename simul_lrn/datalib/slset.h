/***********************************************************************/
/* file slset.h                                                        */
/* contains the class definition of class slset.                       */
/***********************************************************************/

#ifndef SLSET_H
#define SLSET_H

#include "base_slist.h"
#include "stack.h"

/***********************************************************************/
/* class definition of class slset. slset is a set implemented as a    */
/* singli linked (sl) list                                             */
/***********************************************************************/

template <class type>
class slset: public base_slist<type>
{
public:
   slset() { }                                     // constructor
   slset(const slset<type>&) { }                   // copy constructor
   const slset<type>& operator=                    // copy operator
     (const slset<type>&);
   ~slset() { }                                    // destructor
   virtual slist_item<type> *add(const type&);     // prepends
   virtual int card();                             // returns size of slset
   slset<type>& operator|=(const slset<type>&);    // unification
   slset<type>& operator&=(const slset<type>&);    // intersection
   slset<type>& operator-=(const slset<type>&);    // substraction
   int operator==(const slset<type>&);             // test equality of elements
};

#include "slset.cc"                                // because of templates

#endif
