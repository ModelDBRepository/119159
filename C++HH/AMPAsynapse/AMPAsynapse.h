#ifndef _AMPAsynapse_H_
#define _AMPAsynapse_H_
#define A_DIM 1
#define A_PAR 7
class AMPAsynapse:public Synapse{
public:
double t,*x,*dx,*parameter,*extra,*xnew;
int steptransfer;
AMPAsynapse();
AMPAsynapse(Neuron *,Neuron *,double,double,int);
AMPAsynapse(Neuron *, Neuron *,double *a,double *b,int);
~AMPAsynapse();
void setvalues(Neuron *, Neuron *,double *a,double *b,int);
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
