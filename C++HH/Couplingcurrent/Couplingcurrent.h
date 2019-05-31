#ifndef _Couplingcurrent_H_
#define _Couplingcurrent_H_
#define I_DIM 2
#define I_PAR 2
class Couplingcurrent:public Synapse{

public:
double *x,*parameter;
Couplingcurrent();
Couplingcurrent(Neuron *,Neuron *,double );
~Couplingcurrent();
void setvalues(Neuron *,Neuron *,double);
virtual double Isyn();
};
#endif
