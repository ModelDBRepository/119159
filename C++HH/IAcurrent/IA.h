/*Role of IA currents */

/*
IA currents play role in getting continuous gain function for spiking neuron...
Also in there presence there is significant delay between firing of neuron in
response to a depolarizing dc current pulse
*/

#ifndef _IAcurrent_H_
#define _IAcurrent_H_
#define IA_DIM 5
#define IA_PAR 5
class IA:public Ionchannel{

public:
double t,*x,*dx,*xnew,*parameter,*extra;
IA();
IA(Neuron *,double a );
IA(Neuron *,double *a,double *b );
~IA();
void setvalues(Neuron *,double *a,double *b);
void update(double timestep);
 void steptime();
void intime(double intim){t=intim;}
void getvol();
virtual double Iion();
};
#endif
