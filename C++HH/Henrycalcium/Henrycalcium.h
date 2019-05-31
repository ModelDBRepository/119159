#ifndef _Henrycalciumcurrent_H_
#define _Henrycalciumcurrent_H_
#define HC_DIM 1
#define HC_PAR 5
class Henrycalcium:public Calciumchannel{

public:
int Nmda;
double t,*x,*dx,*xnew,*parameter,*extra;
Henrycalcium();
Henrycalcium(Neuron *,Neuron *,double a,int ii);
Henrycalcium(Neuron *,Neuron *,double *a,double *b,int ii );
~Henrycalcium();
void setvalues(Neuron *,Neuron *,double *a,double *b,int ii);
void update(double timestep);
 void steptime();
void intime(double intim){t=intim;}
void getvol();
 double Calcium();
};
#endif
