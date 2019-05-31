/***********************************************************************/
/* file container.h                                                    */
/* contains the class definition of class container.                   */
/***********************************************************************/

#ifndef CONTAINER_H
#define CONTAINER_H

#include <iostream>
#include "container_item.h"

template <class type>
class container;                                   //forward declaration

#include "container_iter.h"

/***********************************************************************/
/* class definition of class container. container is the a generic     */
/* class for all sorts of containers (lists, trees, ...)               */
/***********************************************************************/

template <class type>
class container
{
protected:
   int sz;                                         // contains the size

public:
   container();                                    // constructor
   container(const container&);
   const container<type>& operator=(const container<type>&); // copy operator
   virtual ~container() { }                        // destructor
   virtual container_item<type> *add(const type&)= 0;     // prepends
   virtual int in(const type&)= 0;                        // is arg in the cont. ?
   virtual int del(const type&)= 0;                       // searches & deletes
   int empty();                                    // is container empty ?
   int size();                                     // returns size of container
   virtual cont_iterator<type> *iterator()= 0;     // returns a new iterator
};

#include "container_item.cc"
#include "container.cc"

#endif
