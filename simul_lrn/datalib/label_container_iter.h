/***********************************************************************/
/* file label_container_iter.h                                               */
/* contains the class definition of class label_container_iter which is      */
/* an abstract iterator class for label_container.                           */
/*                                                                     */
/* warning: this class is only a auxiliary class for class label_container.  */
/* NEVER USE THIS CLASS BY ITSELF.                                     */
/* exception: use with another data structure, then make sure that     */
/* label_container_iter.cc is also included.                                 */
/***********************************************************************/

#ifndef LABEL_CONTAINER_ITER_H
#define LABEL_CONTAINER_ITER_H

#include "label_container_item.h"
#include "label_container.h"

#define forall(X) X->init(); while(X->step())

template <class type>
class label_cont_iterator
{
public:
   label_cont_iterator() { }                             // default constr.
   label_cont_iterator(const label_container<type> *) { }       // constructor
   label_cont_iterator(const label_cont_iterator<type> &) { }   // copy constructor
   virtual ~label_cont_iterator() { }                     // destructor
   virtual void init()= 0;                          // initialize
   virtual int step()= 0;                           // step to next node
   virtual label_container_item<type> *current()= 0;      // return thecurr node
   virtual type c_value()= 0;
   virtual int c_label()= 0;
};

#endif
