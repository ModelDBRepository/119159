/***********************************************************************/
/* file base_slist.h                                                   */
/* contains the class definition of class base_slist.                  */
/***********************************************************************/

#ifndef BASE_SLIST_H
#define BASE_SLIST_H

template <class type>
class base_slist;                                // forward declaration

#include "container.h"
#include "slist_item.h"
#include "slist_iter.h"


/***********************************************************************/
/* class definition of class base_slist. base_slist is the basic class */
/* other slist classes are build upon.                                 */
/***********************************************************************/

template <class type>
class base_slist: public container<type> 
{
friend class slist_iterator<type>;               // slist_iterator may access
                                                 // private data fields
protected:
   slist_item<type> *head, *tail;                // points start/end of slist
   slist_item<type> *parent_search(const type&); // internal search routine
                                                 // returns also parent node
public:
   base_slist();                                 // constructor
   base_slist(const base_slist<type>&);          // copy constructor
   const base_slist<type>& operator=             // copy operator
     (const base_slist<type>&);
   virtual ~base_slist();                        // destructor
   virtual slist_item<type> *add(const type&);   // prepends
   virtual int in(const type&);                  // is arg in the list ?
   virtual int del(const type&);                 // searches & deletes
   virtual void clear();                         // empty the list
   virtual slist_iterator<type> *iterator();     // returns a new iterator
};

#include "slist_item.cc"
#include "slist_iter.cc"
#include "base_slist.cc"                         // because of templates
      
#endif
