/****************************************************************************/
/* Class header tri_tnvector.h   (triangle vector)                            */
/* This file is the header file of a 2-dim vector with the property that    */
/* it is automatically antisymmetric, changing the value of (x,y)           */
/* automatically changes (y,x) to -(x,y)!                                   */
/* The type must be a numeric type which supplies the prefix minus.         */
/****************************************************************************/

#ifndef TRI_VECTOR_H
#define TRI_VECTOR_H

template <class type> 
class tri_vector
{
private:
    int size;
    type* vec;

public:
    tri_vector(int);                      // constructor 
    tri_vector(const tri_vector<type>&);  // copy constructor
    const tri_vector<type> &operator=     // copy operator
      (const tri_vector<type> &);  
    ~tri_vector();                        // destructor
    type get(int, int);                   // returns the content of passed pos.
    void set(int, int, type);             // sets a value at passed pos.
};

#include "tri_vector.cc"                  // necessary because of templates

#endif
