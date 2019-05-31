/***********************************************************************/
/* file label_container.h                                                    */
/* contains the class definition of class label_container.                   */
/***********************************************************************/

#ifndef LABEL_CONTAINER_H
#define LABEL_CONTAINER_H

#include <iostream>
#include "label_container_item.h"

template <class type>
class label_container;                                   //forward declaration

#include "label_container_iter.h"

/***********************************************************************/
/* class definition of class label_container. label_container is the a generic     */
/* class for all sorts of label_containers (lists, trees, ...)               */
/***********************************************************************/

template <class type>
class label_container
{
protected:
   int sz;                                         // contains the size

public:
   label_container();                                    // constructor
   label_container(const label_container&);
   const label_container<type>& operator=(const label_container<type>&); // copy operator
   virtual ~label_container() { }                        // destructor
   virtual label_container_item<type> *add(int, const type&)= 0;     // prepends
   virtual int in(const type&)= 0;                        // is arg in the cont. ?
   virtual int del(const type&)= 0;                       // searches & deletes
   int empty();                                    // is label_container empty ?
   int size();                                     // returns size of label_container
   virtual label_cont_iterator<type> *iterator()= 0;     // returns a new iterator
};

#include "label_container_item.cc"
#include "label_container.cc"

#endif
