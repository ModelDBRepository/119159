/**********************************************************************
 file label_container_item.h                                               
 contains the class definition of class label_container item.              
                                                                     
 warning: this class is only a auxiliary class for class label_container.  
 NEVER USE THIS CLASS BY ITSELF.                                     
 exception: use with another data structure, then make sure that     
 label_container_item.cc is also included.                                 
**********************************************************************/

#ifndef LABEL_CONTAINER_ITEM_H
#define LABEL_CONTAINER_ITEM_H

/**********************************************************************
 class definition of class label_container_item. label_container_item
 is the building block for label_container, generic class.       
**********************************************************************/

template <class type>
class label_container_item
{
 public:
   label_container_item(int, const type&);                    // constructor
   label_container_item(const label_container_item &);   // copy constructor
   const label_container_item& operator=                 // copy operator
     (const label_container_item &);
   virtual ~label_container_item() { }                   // destructor
   int label;                                      // contains a label 
   type data;                                      // contains the data
};

#endif
