#ifndef _TwoDsynapse_H_
#define _TwoDsynapse_H_
#define Two_DIM 2
#define Two_PAR 7
class TwoDsynapse:public Synapse{
public:
double t,*x,*dx,*parameter,*extra,*xnew;
int steptransfer;
TwoDsynapse();
TwoDsynapse(Neuron *,Neuron *,double,double,int);
TwoDsynapse(Neuron *, Neuron *,double *a,double *b,int);
~TwoDsynapse();
void setvalues(Neuron *, Neuron *,double *a,double *b,int);
void update(double timestep);
 void steptime();
void intime(double intim){t=intim;}
double Recordtime(double intime){double out=intime;return out;}
double synapse(){return x[0];}
void getvol(double);
virtual double Isyn();
virtual double Calcium();

};
#endif
