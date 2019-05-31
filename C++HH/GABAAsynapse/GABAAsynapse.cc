

#include "iostream.h"
//#include "../define.h"
#include "../Neuron.h"
#include "../Synapse.h"
//#include "../HHneuron/HHneuron.h"
#include "../Izhineuron/Izhineuron.h"
#include "modelGABAAsynapse.h"
#include"../rk4.h"
#include "GABAAsynapse.h"
#ifndef _GABAAsynapse_
#define  _GABAAsynapse_

GABAAsynapse::GABAAsynapse():Synapse(){}
GABAAsynapse::GABAAsynapse(Neuron *obj1,Neuron *obj2,double a,int s):Synapse(obj1,obj2,GA_DIM,GA_PAR)
{
Evalid=Isynvalid=0;

steptransfer=s;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
extra[0]=double(steptransfer);
extra[1]=0;
x[0]=.00;
parameter[0]=a;
parameter[1]=-82;
parameter[2]=.71429;
//parameter[2]=.05;
parameter[3]=.18;
parameter[4]=2;
parameter[5]=5;
}

GABAAsynapse::GABAAsynapse(Neuron *obj1,Neuron *obj2,double *a,double *b,int s):Synapse(obj1,obj2,GA_DIM,GA_PAR)
{
Evalid=Isynvalid=0;
int i;
steptransfer=s;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
extra[0]=double(steptransfer);
extra[1]=0;
for(i=0;i<Dimno;i++)
x[i]=a[i];
for(i=0;i<Parno;i++) //parameter 0 is the synaptic strength and parameter 1 is
parameter[i]=b[i];  // Ereversal for the synapse,parameter 2 is Vth of GABAAsynapse
}                   //parameter 3 is synaptic timescale for synapse  
GABAAsynapse::~GABAAsynapse()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}
void GABAAsynapse::setvalues(Neuron *obj1, Neuron *obj2,double *a,double *b,int s)
{
int i;
steptransfer=s;
Evalid=Isynvalid=0;
source=obj1;
target=obj2;
source->axon.add(this);
target->den.add(this);

Dimno=GA_DIM;
Parno=GA_PAR;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
extra[0]=double(steptransfer);
extra[1]=0;
for(i=0;i<Dimno;i++)
x[i]=a[i];
for(i=0;i<Parno;i++) //parameter 0 is the synaptic strength and parameter 1 is
parameter[i]=b[i];
}
void GABAAsynapse::update(double timestep)
{
step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelGABAAsynapse);
Isynvalid=1;
}
 void GABAAsynapse::steptime()
{
assert(Isynvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=Isynvalid=0;
}
void GABAAsynapse::getvol()
{
if(!Evalid)
{
if (extra[0]==0)
extra[1]=source->E();
else
extra[1]=source->E();
}
Evalid=1;
}
double GABAAsynapse::Isyn()
{
double Isync;
Isync=-parameter[0]*x[0]*(target->E()-parameter[1]);
return(Isync);
}
#endif

/*
x[0]=the synaptic open channel probability
parameter[0]=synaptic max cinductance
parameter[1]=synaptic reversal potential
parameter[2]= alpha term for rise time for open channel
parameter[3]=beta term for fall time of open channel
parameter[4]=Vp..the threshold for activation of open channel probability
(default=2)
parameter[5]=Kp---normalizing factor (default 5)

*/

