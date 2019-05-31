#include "iostream.h"
#include "../Neuron.h"
#include "../Ionchannel.h"
#include "../HHneuron/HHneuron.h"
#include "modelIA.h"
#include"../rk4.h"
#include "IA.h"
IA::IA():Ionchannel(){}
IA::IA(Neuron *obj2,double a):Ionchannel(obj2,IA_DIM,IA_PAR)
{
Evalid=Iionvalid=0;
int i;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
x[0]=.02;
x[1]=.34;

parameter[0]=a;
parameter[1]=-100;
}
IA::IA(Neuron *obj2,double *a,double *b):Ionchannel(obj2,IA_DIM,IA_PAR)
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
IA::~IA()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}
void IA::setvalues(Neuron *obj2,double *a,double *b)
{
int i;
Evalid=Iionvalid=0;
tar=obj2;
tar->ion.add(this);
Dimno=IA_DIM;
Parno=IA_PAR;
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
void IA::update(double timestep)
{

step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelIA);
Iionvalid=1;
}
 void IA::steptime()
{
assert(Iionvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=Iionvalid=0;
}
void IA::getvol()
{
if(!Evalid)
{
extra=new double[2];
extra[1]=tar->E();
}
Evalid=1;
}
double IA::Iion()
{
double Iiona;

Iiona=-parameter[0]*x[0]*x[1]*(tar->E()-parameter[1]);
//Iiona=-parameter[0]*(tar->E()-parameter[1])*(.6*pow(x[0],4.)*x[2]+.4*pow(x[1],4.)*x[3]);
return(Iiona);
}

/*
parameter[0]=max conductance of the channel
parameter[1]=reversal potential of the channel
x[0]=m1 activation variable
x[1]=m2 activation variable
x[2]=h1 inactivation variable
x[3]=h2 inactivation variable
*/
