#ifndef _GABABsynapse_H_
#define _GABABsynapse_H_
#define GB_DIM 2
#define GB_PAR 15
class GABABsynapse:public Synapse{
public:
double t,*x,*dx,*parameter,*extra,*xnew;
GABABsynapse();
GABABsynapse(Neuron *,Neuron *,double);
GABABsynapse(Neuron *, Neuron *,double *a,double *b);
~GABABsynapse();
void setvalues(Neuron *, Neuron *,double *a,double *b);
void update(double timestep);
 void steptime();
void intime(double intim){t=intim;}
double time(){return t;}
double synapse(){return x[0];}
void getvol();
virtual double Isyn();
};
#endif
