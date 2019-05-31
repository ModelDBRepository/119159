#include "iostream.h"
#include "../Neuron.h"
#include "../Ionchannel.h"
#include "../HHneuron/HHneuron.h"
#include "modelIM.h"
#include"../rk4.h"
#include "IM.h"
IM::IM():Ionchannel(){}
IM::IM(Neuron *obj2,double a):Ionchannel(obj2,IM_DIM,IM_PAR)
{
Evalid=Iionvalid=0;

xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
x[0]=.02;

parameter[0]=a;
parameter[1]=-100;
}
IM::IM(Neuron *obj2,double *a,double *b):Ionchannel(obj2,IM_DIM,IM_PAR)
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
for(i=0;i<Parno;i++) 
parameter[i]=b[i];  
}                     
IM::~IM()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}
void IM::setvalues(Neuron *obj2,double *a,double *b)
{
int i;
Evalid=Iionvalid=0;
tar=obj2;
tar->ion.add(this);
Dimno=IM_DIM;
Parno=IM_PAR;
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
void IM::update(double timestep)
{

step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelIM);
Iionvalid=1;
}
 void IM::steptime()
{
assert(Iionvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=Iionvalid=0;
}
void IM::getvol()
{
if(!Evalid)
{
extra=new double[2];
extra[1]=tar->E();
}
Evalid=1;
}
double IM::Iion()
{
double Iiona;
Iiona=-parameter[0]*(tar->E()-parameter[1])*pow(x[0],1);
return(Iiona);
}

/*
parameter[0]=max conductance of the channel
parameter[1]=reversal potentIMl of the channel
x[0]=m1 activation varIMble

*/
