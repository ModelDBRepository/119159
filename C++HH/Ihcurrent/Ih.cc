#include "iostream.h"
#include "../Neuron.h"
#include "../Ionchannel.h"
#include "../HHneuron/HHneuron.h"
#include "modelIh.h"
#include"../rk4.h"
#include "Ih.h"
Ih::Ih():Ionchannel(){}
Ih::Ih(Neuron *obj2,double a):Ionchannel(obj2,Ih_DIM,Ih_PAR)
{
Evalid=Iionvalid=0;
int i;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
x[0]=0.;


parameter[0]=a;
parameter[1]=-43;
}
Ih::Ih(Neuron *obj2,double *a,double *b):Ionchannel(obj2,Ih_DIM,Ih_PAR)
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
Ih::~Ih()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}
void Ih::setvalues(Neuron *obj2,double *a,double *b)
{
int i;
Evalid=Iionvalid=0;
tar=obj2;
tar->ion.add(this);
Dimno=Ih_DIM;
Parno=Ih_PAR;
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
void Ih::update(double timestep)
{

step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelIh);
Iionvalid=1;
}
 void Ih::steptime()
{
assert(Iionvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=Iionvalid=0;
}
void Ih::getvol()
{
if(!Evalid)
{
extra=new double[2];
extra[1]=tar->E();
}
Evalid=1;
}
double Ih::Iion()
{
double Iiona;

Iiona=-parameter[0]*x[0]*(tar->E()-parameter[1]);
//Iiona=-parameter[0]*(tar->E()-parameter[1])*(.6*pow(x[0],4.)*x[2]+.4*pow(x[1],4.)*x[3]);
return(Iiona);
}

/*
parameter[0]=max conductance of the channel
parameter[1]=reversal potential of the channel
x[0]=m1 activation variable; hyperpolarization
activation inward current
*/
