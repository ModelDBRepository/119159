/****************************************************************************/
/* Class header tri_tnvector.h   (triangle vector)                            */
/* This file is the header file of a 2-dim vector with the property that    */
/* it is automatically antisymmetric, changing the value of (x,y)           */
/* automatically changes (y,x) to -(x,y)!                                   */
/* The type must be a numeric type which supplies the prefix minus.         */
/****************************************************************************/

#ifndef COMP_TRI_VECTOR_H
#define COMP_TRI_VECTOR_H

static int pat[4][4]= {{64,16,4,1},{0,0,0,0},{128,32,8,2},{63,207,243,252}};

class comp_tri_vector
{
private:
    int size;
    unsigned char *vec;

public:
    comp_tri_vector(int);     // constructor
    ~comp_tri_vector();       // destructor
    int get(int, int);        // returns the content of passed pos.
    void set(int, int, int);  // sets a value at passed pos.
};

// #include "comp_tri_vector.cc"     // necessary because of templates

#endif
