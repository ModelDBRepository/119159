#ifndef _STFDsynapse_H_
#define _STFDsynapse_H_
#define A_DIM 5
#define A_PAR 7
class STFDsynapse:public Synapse{
public:
double t,*x,*dx,*parameter,*extra,*xnew;
int steptransfer;
STFDsynapse();
STFDsynapse(Neuron *,Neuron *,double,double,int);
STFDsynapse(Neuron *, Neuron *,double *a,double *b,int);
~STFDsynapse();
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
