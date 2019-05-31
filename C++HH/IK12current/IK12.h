#ifndef _IK12current_H_
#define _IK12current_H_
#define IK12_DIM 5
#define IK12_PAR 5
class IK12:public Ionchannel{

public:
double t,*x,*dx,*xnew,*parameter,*extra;
IK12();
IK12(Neuron *,double a );
IK12(Neuron *,double *a,double *b );
~IK12();
void setvalues(Neuron *,double *a,double *b);
void update(double timestep);
 void steptime();
void intime(double intim){t=intim;}
void getvol();
virtual double Iion();
};
#endif
