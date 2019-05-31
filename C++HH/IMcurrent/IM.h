

#ifndef _IMcurrent_H_
#define _IMcurrent_H_
#define IM_DIM 5
#define IM_PAR 5
class IM:public Ionchannel{

public:
double t,*x,*dx,*xnew,*parameter,*extra;
IM();
IM(Neuron *,double a );
IM(Neuron *,double *a,double *b );
~IM();
void setvalues(Neuron *,double *a,double *b);
void update(double timestep);
 void steptime();
void intime(double intim){t=intim;}
void getvol();
virtual double Iion();
};
#endif
