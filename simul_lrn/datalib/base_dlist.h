/***********************************************************************/
/* file base_dlist.h                                                   */
/* contains the class definition of class base_dlist.                  */
/***********************************************************************/

#ifndef BASE_DLIST_H
#define BASE_DLIST_H

template <class type>
class base_dlist;                                // forward declaration

#include "container.h"
#include "dlist_item.h"
#include "dlist_iter.h"


/***********************************************************************/
/* class definition of class base_dlist. base_dlist is the basic class */
/* other dlist classes are build upon. (dlist == doubly linked list)   */
/***********************************************************************/

template <class type>
class base_dlist: public container<type> 
{
friend class dlist_iterator<type>;               // dlist_iterator may access
                                                 // private data fields
 private:
   dlist_item<type> *find(const type&);          // find a data entry
 
 protected:
   dlist_item<type> *head, *tail;                // points start/end of dlist

 public:
   base_dlist();                                 // constructor
   virtual ~base_dlist();                        // destructor
   virtual dlist_item<type> *add(const type&);   // prepends
   virtual int in(const type&);                  // is arg in the list ?
   virtual void del_item(dlist_item<type> *);    // del item from list
   virtual int del(const type&);                 // searches & deletes
   virtual void clear();                         // empty the list
   virtual dlist_iterator<type> *iterator();     // returns a new iterator
};

#include "dlist_item.cc"
#include "dlist_iter.cc"
#include "base_dlist.cc"                         // because of templates
      
#endif
