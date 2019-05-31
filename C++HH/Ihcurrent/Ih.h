/*Role of Ih currents */

/*
Ih currents play role in getting continuous gain function for spiking neuron...
Also in there presence there is significant delay between firing of neuron in
response to a depolarizing dc current pulse
*/

#ifndef _Ihcurrent_H_
#define _Ihcurrent_H_
#define Ih_DIM 5
#define Ih_PAR 5
class Ih:public Ionchannel{

public:
double t,*x,*dx,*xnew,*parameter,*extra;
Ih();
Ih(Neuron *,double a );
Ih(Neuron *,double *a,double *b );
~Ih();
void setvalues(Neuron *,double *a,double *b);
void update(double timestep);
 void steptime();
void intime(double intim){t=intim;}
void getvol();
virtual double Iion();
};
#endif
