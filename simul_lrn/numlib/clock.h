
#ifndef CLOCK_H
#define CLOCK_H

class clock
{
private:
    int start, ende, ord;
    int start_rel, ende_rel;
    int first_time;
    tnvector<int> cnt;

public:
    clock() {}
    clock(int, int, int, int, int);
    ~clock() {}
    void reset(int, int, int, int, int);
    int advance();
    int operator[](int);
    tnvector<int> all();
    tnvector<int> all_up_ordered();
    tnvector<int> all_down_ordered();
};

#include "clock.cc"

#endif
