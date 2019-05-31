#include "iostream.h"
//#include "../define.h"
#include "../Neuron.h"
#include "../Synapse.h"
#include "../HHneuron/HHneuron.h"
#include "modelHenrysynapse.h"
#include"../rk4.h"
#include "Henrysynapse.h"
Henrysynapse::Henrysynapse():Synapse(){}

Henrysynapse::Henrysynapse(Neuron *obj2,double sp,double a,double b,int ii):Synapse((Neuron *)NULL,obj2,H_DIM,H_PAR)
{
Evalid=Isynvalid=0;


Nmda=ii;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
extra[0]=extra[1]=0;
x[0]=0.00;
spike=sp;
if(!Nmda)
{
parameter[0]=a;
parameter[1]=0;
parameter[2]=1.4;
parameter[3]=15.0/14;
parameter[4]=1;
parameter[5]=b;
}
else
{
parameter[0]=a;
parameter[1]=0;
parameter[2]=67.5;
parameter[3]=70./67.5;
parameter[4]=1;
parameter[5]=b;
}
}

Henrysynapse::Henrysynapse(Neuron *obj2,double sp,double *a,double *b,int ii):Synapse((Neuron *)NULL,obj2,H_DIM,H_PAR)
{
Evalid=Isynvalid=0;
int i;

Nmda=ii;
spike=sp;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
extra[0]=extra[1]=0;
for(i=0;i<Dimno;i++)
x[i]=a[i];
for(i=0;i<Parno;i++) //parameter 0 is the synaptic strength and parameter 1 is
parameter[i]=b[i];  // Ereversal for the synapse,parameter 2 is Vth of Henrysynapse
}                   //parameter 3 is synaptic timescale for synapse  
Henrysynapse::~Henrysynapse()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}
void Henrysynapse::setvalues(Neuron *obj2,double sp,double *a,double
*b,int ii)
{
int i;
Evalid=Isynvalid=0;

Nmda=ii;
spike=sp;
target=obj2;

target->den.add(this);

Dimno=H_DIM;
Parno=H_PAR;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
extra[0]=extra[1]=0;
for(i=0;i<Dimno;i++)
x[i]=a[i];
for(i=0;i<Parno;i++) //parameter 0 is the synaptic strength and parameter 1 is
parameter[i]=b[i];
}
void Henrysynapse::update(double timestep)
{
step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelHenrysynapse);
Isynvalid=1;
}
 void Henrysynapse::steptime()
{
assert(Isynvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=Isynvalid=0;
}
void Henrysynapse::getvol()
{
if(!Evalid)
{
extra[1]=spike;
}
Evalid=1;
}
double Henrysynapse::Isyn()
{
double Isync,B;
if(!Nmda)
B=1;
else 
B = 1. / (1 + .288 * parameter[4] * exp (-.062 * target->E()));

Isync=-parameter[0]*x[0]*B*(target->E()-parameter[1]);

return(Isync);
}

double Henrysynapse::Calcium()
{
double Isync,B;
if(!Nmda)
B=1;
else 
B = 1. / (1 + .288 * parameter[4] * exp (-.062 * target->E()));

Isync=-parameter[5]*x[0]*B*(((75./40)*target->E()-65)-parameter[1]);

return(Isync);
}

void Henrysynapse::set_spike(double a)
{
spike=a;
}

/*
parameter[0]=strength of synapse
parameter[1]=reversal potential
parameter[2]=taun/taua
parameter[3]=S1n/S1a
parameter[4]=Mg conc
*/
