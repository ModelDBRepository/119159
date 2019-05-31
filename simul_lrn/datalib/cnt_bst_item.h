/***********************************************************************/
/* file cnt_bst_item.h                                                 */
/* contains the class definition of class cnt_bst_item.                */
/*                                                                     */
/* warning: this class is only a auxiliary class for class cnt_bst.    */
/* NEVER USE THIS CLASS BY ITSELF.                                     */
/* exception: use with another data structure, then make sure that     */
/* cnt_bst_item.cc is also included.                                   */
/***********************************************************************/

#ifndef CNT_BST_ITEM_H
#define CNT_BST_ITEM_H

#include "base_bst.h"

/***********************************************************************/
/* class definition of class cnt_bst_item. bst_item is the building    */
/* for cnt_bst_set.                                                    */
/***********************************************************************/

template <class type>
class cnt_bst_item: public bst_item<type>
{
 public:
    cnt_bst_item(const type&, int cnt= 1);        // constructor
    cnt_bst_item(const cnt_bst_item<type>&);      // copy constructor
    const cnt_bst_item<type>& operator=           // copy operator
      (const cnt_bst_item<type>&);
    virtual ~cnt_bst_item() { }                   // destructor
    int count;                                    // counts multiplicity
};

#endif
