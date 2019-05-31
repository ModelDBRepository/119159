#ifndef _GABAAsynapse_H_
#define _GABAAsynapse_H_
#define GA_DIM 1
#define GA_PAR 6
class GABAAsynapse:public Synapse{
public:
double t,*x,*dx,*parameter,*extra,*xnew;
int steptransfer;
GABAAsynapse();
GABAAsynapse(Neuron *,Neuron *,double,int);
GABAAsynapse(Neuron *, Neuron *,double *a,double *b,int);
~GABAAsynapse();
void setvalues(Neuron *, Neuron *,double *a,double *b,int);
void update(double timestep);
 void steptime();
void intime(double intim){t=intim;}
double time(){return t;}
double synapse(){return x[0];}
void getvol();
virtual double Isyn();

};
#endif
