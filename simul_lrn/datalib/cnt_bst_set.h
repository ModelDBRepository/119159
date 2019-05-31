/***********************************************************************/
/* file cnt_bst_set.h                                                  */
/* contains the class definition of class cnt_bst_set.                 */
/***********************************************************************/

#ifndef CNT_BST_SET_H
#define CNT_BST_SET_H

#include <assert.h>
#include "cnt_bst.h"

/***********************************************************************/
/* class definition of class cnt_bst_set. cnt_bst_set is a set class   */
/* based on a cnt_bst.                                                 */
/***********************************************************************/

template <class type>
class cnt_bst_set: public cnt_bst<type>
{
   friend class cnt_bst_iterator<type>;          // cnt_bst_iterator may access
                                                 // private data fields
public:
   cnt_bst_set() { }                             // constructor
   cnt_bst_set(const cnt_bst_set<type>&) { }     // copy constructor
   const cnt_bst_set<type>& operator=            // copy operator
     (const cnt_bst_set<type>&);
   virtual ~cnt_bst_set() { }                    // destructor
   int card();                                   // returns the size
   type min();                                   // returns minimum
   type max();                                   // returns maximum
   type mean();                                  // returns the (weighted) mean
   cnt_bst_set<type>& operator|=                 // does union
     (const cnt_bst_set<type>&);
   cnt_bst_set<type>& operator&=                 // does intersection
     (const cnt_bst_set<type>&);
   cnt_bst_set<type>& operator-=                 // does substraction
     (const cnt_bst_set<type>&);
   int operator==(const cnt_bst_set<type>&);     // tests equality
};

#include "cnt_bst_set.cc"                        // because of templates

#endif
