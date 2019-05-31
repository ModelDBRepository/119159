#include "iostream.h"
#include "../Neuron.h"
#include "../Calciumchannel.h"
#include "../HHneuron/HHneuron.h"
#include "modelHenrycalcium.h"
#include"../rk4.h"
#include "Henrycalcium.h"
Henrycalcium::Henrycalcium():Calciumchannel(){}
Henrycalcium::Henrycalcium(Neuron *obj1,Neuron *obj2,double a,int ii):Calciumchannel(obj1,obj2,HC_DIM,HC_PAR)
{
Evalid=Calvalid=0;
int i;
Nmda=ii;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
x[0]=.0024;
if(!Nmda)
{
parameter[0]=a;
parameter[1]=0;
parameter[2]=1.4;
parameter[3]=15./14;
parameter[4]=1;
}
else
{
parameter[0]=a;
parameter[1]=0;
parameter[2]=69.75;
parameter[3]=70./69.75;
parameter[4]=1;
}
}

Henrycalcium::Henrycalcium(Neuron *obj1,Neuron *obj2,double *a,double *b,int ii):Calciumchannel(obj1,obj2,HC_DIM,HC_PAR)
{
Evalid=Calvalid=0;
int i;
Nmda=ii;
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
Henrycalcium::~Henrycalcium()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}
void Henrycalcium::setvalues(Neuron *obj1,Neuron *obj2,double *a,double *b,int ii)
{
int i;
Nmda=ii;
Evalid=Calvalid=0;
sourcecal=obj1;
tarcal=obj2;
sourcecal->calsource.add(this);
tarcal->cal.add(this);
Dimno=HC_DIM;
Parno=HC_PAR;
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
void Henrycalcium::update(double timestep)
{

step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelHenrycalcium);
Calvalid=1;
}
 void Henrycalcium::steptime()
{
assert(Calvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=Calvalid=0;
}
void Henrycalcium::getvol()
{
if(!Evalid)
{
extra=new double[2];
extra[1]=sourcecal->E();
}
Evalid=1;
}
double Henrycalcium::Calcium()
{
double Icala,B;
if(!Nmda)
B=1;
else 
B = 1. / (1 + .288 * parameter[4] * exp (-.062 * tarcal->E()));
Icala=-parameter[0]*x[0]*B*(tarcal->E()-parameter[1]);
return(Icala);
}
/*
parameter[0]=strength of calcium entry through synapse
parameter[1]=reversal potential
parameter[2]=taun/taua
parameter[3]=S1n/S1a
parameter[4]=Mg conc
*/
