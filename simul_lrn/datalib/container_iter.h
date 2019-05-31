/***********************************************************************/
/* file container_iter.h                                               */
/* contains the class definition of class container_iter which is      */
/* an abstract iterator class for container.                           */
/*                                                                     */
/* warning: this class is only a auxiliary class for class container.  */
/* NEVER USE THIS CLASS BY ITSELF.                                     */
/* exception: use with another data structure, then make sure that     */
/* container_iter.cc is also included.                                 */
/***********************************************************************/

#ifndef CONTAINER_ITER_H
#define CONTAINER_ITER_H

#include "container_item.h"
#include "container.h"

#define forall(X) X->init(); while(X->step())

template <class type>
class cont_iterator
{
public:
   cont_iterator() { }                             // default constr.
   cont_iterator(const container<type> *) { }       // constructor
   cont_iterator(const cont_iterator<type> &) { }   // copy constructor
   virtual ~cont_iterator() { }                     // destructor
   virtual void init()= 0;                          // initialize
   virtual int step()= 0;                           // step to next node
   virtual container_item<type> *current()= 0;      // return thecurr node
   virtual type c_value()= 0;
};

#endif
