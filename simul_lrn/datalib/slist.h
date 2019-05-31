/***********************************************************************/
/* file slist.h                                                        */
/* contains the class definition of class slist                        */
/***********************************************************************/

#ifndef SLIST_H
#define SLIST_H

#include <assert.h>
#include "base_slist.h"
#include "slist.cc"                                // because of templates
/***********************************************************************/
/* class definition of class slist. slist is a singly linked list.     */
/***********************************************************************/

template <class type>
class slist: public base_slist<type>
{
   friend class slist_iterator<type>;              // slist_iterator may access
                                                   // private data fields

public:
   slist() { }                                     // constructor
   slist(const slist<type>&) { }                   // copy constructor
   const slist<type>& operator=                    // copy operator
     (const slist<type>&);
   virtual ~slist() { }                            // destructor (empty)
   virtual slist_item<type> *prepend(const type&); // prepends item with add
   virtual type del_item(slist_item<type> *);      // removes arg. from slist
   virtual slist_item<type> *first();              // returns first node
   virtual slist_item<type> *last();               // returns last node
   virtual slist_item<type> *search(const type&);  // searches node 
   virtual slist_item<type> *succ                  // returns successor 
     (slist_item<type> *);
   virtual slist_item<type> *insert                // inserts a node after a
     (slist_item<type> *, const type&);            // specific node 
   virtual slist_item<type> *append(const type&);  // appends node
};

template <class type>
ostream &operator<<(ostream&, slist<type>&);
template <class type>
istream &operator>>(istream&, slist<type>&);



#endif
