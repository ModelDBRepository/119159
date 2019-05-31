#include "iostream.h"
#include "../Neuron.h"
#include "../Calciumchannel.h"
#include "../HHneuron/HHneuron.h"
#include "modelIL.h"
#include"../rk4.h"
#include "IL.h"
IL::IL():Calciumchannel(){}
IL::IL(Neuron *obj2,double a,double b):Calciumchannel((Neuron *) NULL,obj2,IL_DIM,IL_PAR)
{
Evalid=Icalvalid=0,Calvalid=0;
int i;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
x[0]=.00;
x[1]=0.;
parameter[0]=a;
parameter[1]=b;
parameter[2]=100.0;

}
IL::IL(Neuron *obj2,double *a,double *b):Calciumchannel((Neuron *) NULL,obj2,IL_DIM,IL_PAR)
{
Evalid=Icalvalid=0,Calvalid=0;
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
IL::~IL()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}
void IL::setvalues(Neuron *obj2,double *a,double *b)
{
int i;
Evalid=Icalvalid=0,Calvalid=0;
tarcal=obj2;
tarcal->cal.add(this);
Dimno=IL_DIM;
Parno=IL_PAR;
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
void IL::update(double timestep)
{

step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelIL);
Icalvalid=1;
Calvalid=1;
}
 void IL::steptime()
{
assert(Icalvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=Icalvalid=Calvalid=0;
}
void IL::getvol()
{
if(!Evalid)
{
extra=new double[2];
extra[1]=tarcal->E();
}
Evalid=1;
}
double IL::Icalcium()
{
double v,minf,Ical;
v=tarcal->E();
minf=1./(1+exp(-(v+20.0)/2.5));
Ical=-parameter[0]*pow(minf,1.)*(tarcal->E()-parameter[2]);
return(Ical);
}
double IL::Calcium()
{
double v, minf,Ccal;
v=tarcal->E();
minf=1./(1+exp(-(v+25.0)/2.5));
Ccal=-parameter[1]*minf*(tarcal->E()-parameter[2]);
return(Ccal);
}

/*
x[0]= m activation variable
x[1]=h inactivation avriable
parameter[0]= conductance for calcium current in neuron
parameter[1]=Calcium conc conductance term in calcium dynamics
*/
