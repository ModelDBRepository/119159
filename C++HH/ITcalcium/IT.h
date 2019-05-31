#ifndef _ITcurrent_H_
#define _ITcurrent_H_
#define IT_DIM 5
#define IT_PAR 5
class IT:public Calciumchannel{

public:
double t,*x,*dx,*xnew,*parameter,*extra;
int Val;
IT();
IT(Neuron *,double a,double b,int GHK);
IT(Neuron *,double *a,double *b );
~IT();
void setvalues(Neuron *,double *a,double *b);
void update(double timestep);
 void steptime();
void intime(double intim){t=intim;}
void getvol();
 double Icalcium();
double Calcium();
};
#endif
