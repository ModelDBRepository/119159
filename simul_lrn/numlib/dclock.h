
#ifndef DCLOCK_H
#define DCLOCK_H

#include "../lib2/bst_set.h"

class dclock
{
private:
    int start, ende, ord;
    int start_rel, ende_rel;
    int first_time;
    tnvector<int> cnt;

public:
    dclock() {}
    dclock(int, int, int, int, int);
    ~dclock() {}
    void reset(int, int, int, int, int);
    int advance();
    int operator[](int);
    tnvector<int> all();
    tnvector<int> all_up_ordered();
    tnvector<int> all_down_ordered();
};

#include "dclock.cc"

#endif
