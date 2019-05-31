/***********************************************************************/
/* file container_item.h                                               */
/* contains the class definition of class container item.              */
/*                                                                     */
/* warning: this class is only a auxiliary class for class container.  */
/* NEVER USE THIS CLASS BY ITSELF.                                     */
/* exception: use with another data structure, then make sure that     */
/* container_item.cc is also included.                                 */
/***********************************************************************/

#ifndef CONTAINER_ITEM_H
#define CONTAINER_ITEM_H

/***********************************************************************/
/* class definition of class container_item. container_item is the     */
/* building block for container, generic class.                        */
/***********************************************************************/

template <class type>
class container_item
{
 public:
   container_item(const type&);                    // constructor
   container_item(const container_item &);         // copy constructor
   const container_item& operator=                 // copy operator
     (const container_item &);
   virtual ~container_item() { }                   // destructor
   type data;                                      // contains the data
};

#endif
