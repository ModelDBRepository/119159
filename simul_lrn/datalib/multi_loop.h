
#ifndef MULTI_LOOP_H
#define MULTI_LOOP_H

template<class type> class multi_loop
{
private:
  type min, max;
  type *cnt;
  short d, dpt;

public:
  multi_loop();
  multi_loop(type, type, short);
  ~multi_loop();
  void reset();
  int step();
  type operator[](short);
  type *vals();
  short depth();
  void init(type, type, short);
};

#include "multi_loop.cc"

#endif


