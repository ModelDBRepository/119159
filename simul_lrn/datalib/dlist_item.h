/***********************************************************************/
/* file dlist_item.h                                                   */
/* contains the class definition of class dlist_item.                  */
/*                                                                     */
/* warning: this class is only a auxiliary class for class dlist.      */
/* NEVER USE THIS CLASS BY ITSELF.                                     */
/* exception: use with another data structure, then make sure that     */
/* dlist_item.cc is also included.                                     */
/***********************************************************************/


#ifndef DLIST_ITEM_H
#define DLIST_ITEM_H

/***********************************************************************/
/* class definition of class dlist_item. dlist_item is the             */
/* building block for all types of dlists (singly linked lists).       */
/***********************************************************************/

template <class type>
class dlist_item: public container_item<type>
{
 public:
   dlist_item(const type&);                         // constructor
   ~dlist_item();                                   // destructor
   dlist_item<type> *prev;                          // points to predecessor   
   dlist_item<type> *next;                          // points to successor
};


#endif
