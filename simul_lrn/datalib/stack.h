/***********************************************************************/
/* file stack.h                                                        */
/* contains the class definition of class stack.                       */
/***********************************************************************/

#ifndef STACK_H
#define STACK_H

#include <assert.h>
#include "base_slist.h"

/***********************************************************************/
/* class definition of class stack.                                    */
/***********************************************************************/

template <class type>
class stack:public base_slist<type>
{
public:
   stack() { }                                     // constructor
   stack(const stack<type>&) { }                   // copy constructor
   const stack<type>& operator=                    // copy operator
     (const stack<type>&);
   virtual ~stack() { }                            // destructor
   void push(const type&);                         // prepends item with add
   type pop();                                     // removes first from slist
   type peek();                                    // returns value of first
};

#include "stack.cc"                                // because of templates

#endif
