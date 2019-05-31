#ifndef _NMDAsynapse_H_
#define _NMDAsynapse_H_
#define N_DIM 1
#define N_PAR 8
class NMDAsynapse:public Synapse{
public:
double t,*x,*dx,*parameter,*extra,*xnew;
NMDAsynapse();
NMDAsynapse(Neuron *,Neuron *,double,double);
NMDAsynapse(Neuron *, Neuron *,double *a,double *b);
~NMDAsynapse();
void setvalues(Neuron *, Neuron *,double *a,double *b);
void update(double timestep);
 void steptime();
void intime(double intim){t=intim;}
double time(){return t;}
double synapse(){return x[0];}
void getvol();
virtual double Isyn();
virtual double Calcium();
};
#endif
