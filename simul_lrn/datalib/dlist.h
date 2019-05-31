/***********************************************************************/
/* file dlist.h                                                        */
/* contains the class definition of class dlist                        */
/***********************************************************************/

#ifndef DLIST_H
#define DLIST_H

#include <assert.h>
#include "base_dlist.h"

/***********************************************************************/
/* class definition of class dlist. dlist is a singly linked list.     */
/***********************************************************************/

template <class type>
class dlist: public base_dlist<type>
{
   friend class dlist_iterator<type>;              // dlist_iterator may access
                                                   // private data fields

public:
   dlist() { }                                     // constructor
   virtual ~dlist() { }                            // destructor (empty)
   virtual dlist_item<type> *prepend(const type&); // prepends item with add
   virtual dlist_item<type> *first();              // returns first node
   virtual dlist_item<type> *last();               // returns last node
   virtual dlist_item<type> *succ                  // returns successor 
     (dlist_item<type> *);
   dlist_item<type> *search(const type&);          // find a data entry
   virtual dlist_item<type> *insert_after          // inserts a node after a
     (dlist_item<type> *, const type&);            // specific node 
   virtual dlist_item<type> *append(const type&);  // appends node
};

template <class type>
ostream &operator<<(ostream&, dlist<type>&);

template <class type>
istream &operator>>(istream&, dlist<type>&);

#include "dlist.cc"                                // because of templates

#endif
