#include "iostream.h"
//#include "../define.h"
#include "../Neuron.h"
#include "../Synapse.h"
#include "../HHneuron/HHneuron.h"
#include "modelGABABsynapse.h"
#include"../rk4.h"
#include "GABABsynapse.h"
GABABsynapse::GABABsynapse():Synapse(){}
GABABsynapse::GABABsynapse(Neuron *obj1,Neuron *obj2,double a):Synapse(obj1,obj2,GB_DIM,GB_PAR)
{
Evalid=Isynvalid=0;
//int i;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
extra[0]=extra[1]=0;
x[0]=.005;
x[1]=.002;
parameter[0]=a;
parameter[1]=-80;
parameter[2]=.09;
parameter[3]=.0012;
parameter[4]=2;
parameter[5]=5;
parameter[6]=.18;
parameter[7]=.034;
}

GABABsynapse::GABABsynapse(Neuron *obj1,Neuron *obj2,double *a,double *b):Synapse(obj1,obj2,GB_DIM,GB_PAR)
{
Evalid=Isynvalid=0;
int i;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
extra[0]=extra[1]=0;
for(i=0;i<Dimno;i++)
x[i]=a[i];
for(i=0;i<Parno;i++) //parameter 0 is the synaptic strength and parameter 1 is
parameter[i]=b[i];  // Ereversal for the synapse,parameter 2 is Vth of GABABsynapse
}                   //parameter 3 is synaptic timescale for synapse  
GABABsynapse::~GABABsynapse()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}
void GABABsynapse::setvalues(Neuron *obj1, Neuron *obj2,double *a,double *b)
{
int i;
Evalid=Isynvalid=0;
source=obj1;
target=obj2;
source->axon.add(this);
target->den.add(this);

Dimno=GB_DIM;
Parno=GB_PAR;
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
void GABABsynapse::update(double timestep)
{
step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelGABABsynapse);
Isynvalid=1;
}
 void GABABsynapse::steptime()
{
assert(Isynvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=Isynvalid=0;
}
void GABABsynapse::getvol()
{
if(!Evalid)
{
extra[1]=source->E();
}
Evalid=1;
}
double GABABsynapse::Isyn()
{
double Isync;
Isync=-parameter[0]*pow(x[1],4.)/(5.+pow(x[1],4.))*(target->E()-parameter[1]);
return(Isync);
}
/*
x[0]=the synaptic open channel probability
parameter[0]=synaptic max cinductance
parameter[1]=synaptic reversal potential
parameter[2]= alpha term for rise time for open channel
parameter[3]=beta term for fall time of open channel
parameter[4]=Vp..the threshold for activation of open channel probability
(default=2)
parameter[5]=Kp---normalizing factor (default 5)
parameter[6]=K3 term for x[1] dynamics
parameter[7]=K4 term for x[1] dynamics
*/

