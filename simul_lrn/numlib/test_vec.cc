#include <iostream>

#include "tnvector.h"

int main(void)
{
  tnvector<int> a(3);

  cin >> a;
  for (int i=1; i< 3; i++)
  {
    a= i*(a+a);
    cout << a;
    cout << endl;
  }
  return 0;
}
