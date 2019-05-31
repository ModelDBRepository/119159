/*Role of IA currents */

/*
IA currents play role in getting continuous gain function for spiking neuron...
Also in there presence there is significant delay between firing of neuron in
response to a depolarizing dc current pulse
*/

#ifndef _PDrule_H_
#define _PD_H_
#define PD_DIM 6
#define PD_PAR 8
class PD{

public:
double t,*x,*dx,*xnew,*parameter,*extra,gvalid;
PD();
PD(double *a, double *b);
~PD();
void update(double timestep);
 void steptime();
void intime(double intim){t=intim;}
void getcal(double &,double &);
double deltag(){return x[2];}
void getdiff(double &);
void operator = (PD &test);
};
#endif
