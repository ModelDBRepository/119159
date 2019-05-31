#include <iostream>

#include "matrix.h"

int main(void)
{
  matrix<double> a(3,3);
  matrix<double> c(3,3);
  
  cin >> a;
  cout << a;
  cin >> c;
  cout << c;
  cout << endl;
  cout << a*c;
  cout << endl << endl;


  matrix<double> b(3,3);
  if (a.invert(b))
  {
    cout << endl;
    cout << endl;
    cout << b;
    cout << a;
    cout << a*b;
  }
  else
  {
    cout << "no inverse";
  }
}
