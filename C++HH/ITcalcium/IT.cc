#include "iostream.h"
#include "../Neuron.h"
#include "../Calciumchannel.h"
#include "../HHneuron/HHneuron.h"
#include "modelIT.h"
#include"../rk4.h"
#include "IT.h"
IT::IT():Calciumchannel(){}
IT::IT(Neuron *obj2,double a,double b,int GHK):Calciumchannel((Neuron *) NULL,obj2,IT_DIM,IT_PAR)
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
Val=GHK;
}
IT::IT(Neuron *obj2,double *a,double *b):Calciumchannel((Neuron *) NULL,obj2,IT_DIM,IT_PAR)
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
IT::~IT()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}
void IT::setvalues(Neuron *obj2,double *a,double *b)
{
int i;
Evalid=Icalvalid=0,Calvalid=0;
tarcal=obj2;
tarcal->cal.add(this);
Dimno=IT_DIM;
Parno=IT_PAR;
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
void IT::update(double timestep)
{

step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelIT);
Icalvalid=1;
Calvalid=1;
}
 void IT::steptime()
{
assert(Icalvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=Icalvalid=Calvalid=0;
}
void IT::getvol()
{
if(!Evalid)
{
extra=new double[2];
extra[1]=tarcal->E();
}
Evalid=1;
}
double IT::Icalcium()
{
double Icala,Caout,Co,temp,tabs,faraday,rgas,GHK,xi;
Caout=40000;
Co=1;
temp=25;
tabs=273.15;
faraday=96485;
rgas=8.3147;
xi=tarcal->E()*faraday*2/(rgas*(tabs+temp)*1000);
if (Val)
//GHK=0;
GHK=(tarcal->E()/Co)*(tarcal->Cain()-Caout*exp(-xi))/(1-exp(-xi));
else 
GHK=(tarcal->E()/Co)*(1.-Caout*exp(-xi))/(1-exp(-xi));
Icala=-parameter[0]*GHK*pow(x[0],2.)*x[1];
return(Icala);
}
double IT::Calcium()
{
double Icala,Caout,Co,temp,tabs,faraday,rgas,GHK,xi;
Caout=15000;
Co=1;
temp=25;
tabs=273.15;
faraday=96485;
rgas=8.3147;
xi=tarcal->E()*faraday*2/(rgas*(tabs+temp)*1000);
if (Val)
//GHK=0;
GHK=(tarcal->E()/Co)*(tarcal->Cain()-Caout*exp(-xi))/(1-exp(-xi)); 
//currently some problem with tarcal-Cain() some how it is giving me an error

else
GHK=(tarcal->E()/Co)*(1-Caout*exp(-xi))/(1-exp(-xi));
Icala=-parameter[1]*GHK*pow(x[0],2.)*x[1];
return(Icala);
}

/*
x[0]= m activation variable
x[1]=h inactivation avriable
parameter[0]= conductance for calcium current in neuron
parameter[1]=Calcium conc conductance term in calcium dynamics
*/
