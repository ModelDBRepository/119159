/***********************************************************************/
/* file queue.h                                                        */
/* contains the class definition of class queue                        */
/***********************************************************************/

#ifndef QUEUE_H
#define QUEUE_H

#include <assert.h>
#include "base_slist.h"

/***********************************************************************/
/* class definition of class queue. queue is a singly linked list.     */
/***********************************************************************/

template <class type>
class queue: public base_slist<type>
{
   friend class slist_iterator<type>;              // slist_iterator may access
                                                   // private data fields
public:
   queue() { }                                     // constructor
   queue(const queue<type>&) { }                   // copy constructor
   const queue<type>& operator=                    // copy operator
     (const queue<type>&);
   virtual ~queue() { }                            // destructor (empty)
   void push(const type&);                         // prepends item with add
   type pop();                                     // returns last nd & del. it
   type peek();                                    // returns "last" node
   type last();                                    // returns latest added node
};

#include "queue.cc"                                // because of templates

#endif



