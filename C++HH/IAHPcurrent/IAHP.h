#ifndef _IAHPcurrent_H_
#define _IAHPcurrent_H_
#define IAHP_DIM 5
#define IAHP_PAR 5
class IAHP:public Calciumchannel{

public:
double t,*x,*dx,*xnew,*parameter,*extra;
int Val;
IAHP();
IAHP(Neuron *,double a,double b,double c,double d);
IAHP(Neuron *,double *a);
~IAHP();
void setvalues(Neuron *,double *a);
void update(double timestep);
 void steptime();
void intime(double intim){t=intim;}
void getvol();
 double Icalcium();
 double Calcium();
};
#endif
