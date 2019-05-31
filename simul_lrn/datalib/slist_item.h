/***********************************************************************/
/* file slist_item.h                                                   */
/* contains the class definition of class slist_item.                  */
/*                                                                     */
/* warning: this class is only a auxiliary class for class slist.      */
/* NEVER USE THIS CLASS BY ITSELF.                                     */
/* exception: use with another data structure, then make sure that     */
/* slist_item.cc is also included.                                     */
/***********************************************************************/


#ifndef SLIST_ITEM_H
#define SLIST_ITEM_H

/***********************************************************************/
/* class definition of class slist_item. slist_item is the             */
/* building block for all types of slists (singly linked lists).       */
/***********************************************************************/

template <class type>
class slist_item: public container_item<type>
{
 public:
   slist_item(const type&);                         // constructor
   slist_item(const slist_item<type>&);             // copy constructor
   const slist_item<type>& operator=                // copy operator
     (const slist_item<type>&);
   ~slist_item();                                   // destructor
   slist_item<type> *next;                          // points to successor
};


#endif
