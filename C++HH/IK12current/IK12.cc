#include "iostream.h"
#include "../Neuron.h"
#include "../Ionchannel.h"
#include "../HHneuron/HHneuron.h"
#include "modelIK12.h"
#include"../rk4.h"
#include "IK12.h"
IK12::IK12():Ionchannel(){}
IK12::IK12(Neuron *obj2,double a):Ionchannel(obj2,IK12_DIM,IK12_PAR)
{
Evalid=Iionvalid=0;
int i;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
x[0]=.0033;
x[1]=.23;
x[2]=.21;
parameter[0]=a;
parameter[1]=-95;

}
IK12::IK12(Neuron *obj2,double *a,double *b):Ionchannel(obj2,IK12_DIM,IK12_PAR)
{
Evalid=Iionvalid=0;
int i;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];

for(i=0;i<Dimno;i++)
x[i]=a[i];
for(i=0;i<Parno;i++) //parameter 0 is the synaptic strength and parameter 1 is
parameter[i]=b[i];  // Ereversal for the synapse,parameter 2 is Vth of Rallsynapse
}                   //parameter 3 is synaptic timescale for synapse  
IK12::~IK12()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}
void IK12::setvalues(Neuron *obj2,double *a,double *b)
{
int i;
Evalid=Iionvalid=0;
tar=obj2;
tar->ion.add(this);
Dimno=IK12_DIM;
Parno=IK12_PAR;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
for(i=0;i<Dimno;i++)
x[i]=a[i];
for(i=0;i<Parno;i++) //parameter 0 is the synaptic strength and parameter 1 is
parameter[i]=b[i];
}
void IK12::update(double timestep)
{

step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelIK12);
Iionvalid=1;
}
 void IK12::steptime()
{
assert(Iionvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=Iionvalid=0;
}
void IK12::getvol()
{
if(!Evalid)
{
extra=new double[2];
extra[1]=tar->E();
}
Evalid=1;
}
double IK12::Iion()
{
double Iiona;
Iiona=-parameter[0]*(tar->E()-parameter[1])*x[0]*(.4*x[1]+.6*x[2]);
return(Iiona);
}

/*
parameter[0]=max conductance of the channel
parameter[1]=reversal potential of the channel
x[0]=activation variable
x[1]=1st inactivation variable
x[2]=2nd inactivation variable
*/
