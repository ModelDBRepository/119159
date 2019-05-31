#ifndef _Henrysynapse_H_
#define _Henrysynapse_H_
#define H_DIM 1
#define H_PAR 6
class Henrysynapse:public Synapse{
public:
int Nmda;
double t,*x,*dx,*parameter,*extra,*xnew,spike;
Henrysynapse();
Henrysynapse(Neuron *,double ,double ,double,int);
Henrysynapse(Neuron *,double  ,double *a,double *b,int);
~Henrysynapse();
void setvalues(Neuron *,double ,double *a,double *b,int);
void update(double timestep);
 void steptime();
void intime(double intim){t=intim;}
double time(){return t;}
double synapse(){return x[0];}
void getvol();
void set_spike(double a);
double get_spike(){return spike;}
virtual double Isyn();
virtual double Calcium();
};
#endif
