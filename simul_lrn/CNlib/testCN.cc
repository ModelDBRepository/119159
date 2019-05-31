/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

using namespace std;

#include <iostream>

#include "CN_Poissonneuron.h"
#include "CN_ECneuron.h"
#include "CN_absynapse.h"
#include "CN_NeuronModel.h"
#include "CN_rk65n.h"
#include "slist.h"
#include "CN_DCInput.h"

#include "CN_Poissonneuron.cc"
#include "CN_ECneuron.cc"
#include "CN_absynapse.cc"
#include "CN_NeuronModel.cc"
#include "CN_DCInput.cc"

int main(int argc, char *argv[])
{
  slist<neuron *> neurs;
  slist<synapse *> syns;
  neuron *n, *n2;
  synapse *s, *ins;
  double *x, *xn, *tmp;
  int N;

  double IDC= atof(argv[1]);
  tnvector<int> pos(3,"0 0 0");
  n= new Poissonneuron(0,pos);
  //neurs.append(n);
  n2= new ECneuron(1,pos);
  neurs.append(n2);
  s= new absynapse(n, n2, 0.5, 0.0, -20.0, 0.1, 0.2, 5.0);
  syns.append(s);
  //  ins= new DCInput(n, IDC);
  //  syns.append(ins);
  
  NeuronModel model(&neurs, &syns, N, cerr);
  x= new double[N];
  xn= new double[N];
  //n->init(x, ECN_INIVARS);
  n2->init(x, ECN_INIVARS);
  s->init(x, ABSYN_INIVARS);
  rk65n machine(N, 0.1, 1e-8, 1e-12, 1e-6);

  double dt= 0.1;
  double dtx= 0.1;
  for (int i= 0; i < 100000; i++) {
    dtx= machine.integrate(x, xn, &model, dt);
    cerr << dtx << endl;
    ((Poissonneuron *)n)->validate_E(x, dt);
    dt= dtx;
    tmp= x; x= xn; xn= tmp;
    ((Poissonneuron *)n)->step();
    cout << x[0] << " " << n->E(x) << " " << n2->E(x) << endl;
  }
  return 0;
}
    
