/*********************************************************************
 file base_label_dlist.h                                                   
 contains the class definition of class base_label_dlist.                  
*********************************************************************/

#ifndef BASE_LABEL_DLIST_H
#define BASE_LABEL_DLIST_H

template <class type>
class base_label_dlist;                                // forward declaration

#include "label_container.h"
#include "label_dlist_item.h"
#include "label_dlist_iter.h"


/*********************************************************************
 class definition of class base_label_dlist. base_label_dlist is the
 basic class other label_dlist classes are build upon.
 (label_dlist == doubly linked list)    
*********************************************************************/

template <class type>
class base_label_dlist: public label_container<type> 
{
  friend class label_dlist_iterator<type>; // label_dlist_iterator may access
                                           // private data fields
  friend class us_array<type>;
 
 private:
 
 protected:
   label_dlist_item<type> *head, *tail;     // points start/end of label_dlist

 public:
   base_label_dlist();                                 // constructor
   virtual ~base_label_dlist();                        // destructor
   virtual label_dlist_item<type> *insert
     (label_dlist_item<type> *, int, const type&);
   virtual label_dlist_item<type> *prepend(int, const type &);
   virtual label_dlist_item<type> *add(int, const type&);   // adds
   virtual int in(const type&);                  // is arg in the list ?
   virtual int label_in(int);                  // is arg in the list ?
   virtual void del_item(label_dlist_item<type> *);    // delete passed item
   virtual int del(const type&);                 // searches & deletes
   virtual int label_del(int);                   // searches by label & deletes
   virtual void clear();                         // empty the list
   label_dlist_item<type> *find(const type&); 
   label_dlist_item<type> *label_find(int); 
   type label_find_data(int); 
   virtual label_dlist_iterator<type> *iterator();  // returns a new iterator
};

#include "label_dlist_item.cc"
#include "label_dlist_iter.cc"
#include "base_label_dlist.cc"                         // because of templates
      
#endif
