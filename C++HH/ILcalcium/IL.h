#ifndef _ILcurrent_H_
#define _ILcurrent_H_
#define IL_DIM 5
#define IL_PAR 5
class IL:public Calciumchannel{

public:
double t,*x,*dx,*xnew,*parameter,*extra;
int Val;
IL();
IL(Neuron *,double a,double b);
IL(Neuron *,double *a,double *b);
~IL();
void setvalues(Neuron *,double *a,double *b);
void update(double timestep);
 void steptime();
void intime(double intim){t=intim;}
void getvol();
 double Icalcium();
double Calcium();
};
#endif
