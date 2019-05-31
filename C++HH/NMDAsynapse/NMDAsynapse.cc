#include "iostream.h"
//#include "../define.h"
#include "../Neuron.h"
#include "../Synapse.h"
#include "../HHneuron/HHneuron.h"
#include "modelNMDAsynapse.h"
#include"../rk4.h"
#include "NMDAsynapse.h"
NMDAsynapse::NMDAsynapse():Synapse(){}
NMDAsynapse::NMDAsynapse(Neuron *obj1,Neuron *obj2,double a,double b):Synapse(obj1,obj2,N_DIM,N_PAR)
{
Evalid=Isynvalid=0;

xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
extra[0]=extra[1]=0;
x[0]=.005;
parameter[0]=a;
parameter[1]=0;
parameter[2]=.385715;  //.072;
parameter[3]=.014285;//.0066;
parameter[4]=2;
parameter[5]=5;
parameter[6]=1;
parameter[7]=b;
}

NMDAsynapse::NMDAsynapse(Neuron *obj1,Neuron *obj2,double *a,double *b):Synapse(obj1,obj2,N_DIM,N_PAR)
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
parameter[i]=b[i];  // Ereversal for the synapse,parameter 2 is Vth of NMDAsynapse
}                   //parameter 3 is synaptic timescale for synapse  
NMDAsynapse::~NMDAsynapse()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}
void NMDAsynapse::setvalues(Neuron *obj1, Neuron *obj2,double *a,double *b)
{
int i;
Evalid=Isynvalid=0;
source=obj1;
target=obj2;
source->axon.add(this);
target->den.add(this);

Dimno=N_DIM;
Parno=N_PAR;
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
void NMDAsynapse::update(double timestep)
{
step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelNMDAsynapse);
Isynvalid=1;
}
 void NMDAsynapse::steptime()
{
assert(Isynvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=Isynvalid=0;
}
void NMDAsynapse::getvol()
{
if(!Evalid)
{
extra=new double[2];
extra[1]=source->E();
}
Evalid=1;
}
double NMDAsynapse::Isyn()
{
double Isync,B;
 B = 1. / (1 + .288 * parameter[6] * exp (-.062 * target->E()));
Isync=-parameter[0]*x[0]*B*(target->E()-parameter[1]);
return(Isync);
}
double NMDAsynapse::Calcium()
{
double Isync,B;
 B = 1. / (1 + .288 * parameter[6] * exp (-.062 * target->E()));
Isync=-parameter[7]*x[0]*B*(target->E()-parameter[1]);
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
parameter[6]=Mg conc......(default 1)
*/

