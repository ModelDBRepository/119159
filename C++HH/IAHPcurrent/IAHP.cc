#include "iostream.h"
#include "../Neuron.h"
#include "../Calciumchannel.h"
#include "../HHneuron/HHneuron.h"
#include "modelIAHP.h"
#include"../rk4.h"
#include "IAHP.h"
IAHP::IAHP():Calciumchannel(){}
IAHP::IAHP(Neuron *obj2,double a,double b,double c,double d):Calciumchannel((Neuron *) NULL,obj2,IAHP_DIM,IAHP_PAR)
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
parameter[0]=a; //strength of calcium current conductance
parameter[1]=b;// reversal potential for potassium
parameter[2]=c;// strenght of conductance to model calcium dynamics
parameter[3]=d;// reversal potential to model calcium dynamics

}
IAHP::IAHP(Neuron *obj2,double *a):Calciumchannel((Neuron *) NULL,obj2,IAHP_DIM,IAHP_PAR)
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
}                  
IAHP::~IAHP()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}
void IAHP::setvalues(Neuron *obj2,double *a)
{
int i;
Evalid=Icalvalid=0,Calvalid=0;
tarcal=obj2;
tarcal->cal.add(this);
Dimno=IAHP_DIM;
Parno=IAHP_PAR;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
for(i=0;i<Dimno;i++)
x[i]=a[i];
}
void IAHP::update(double timestep)
{

step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelIAHP);
Icalvalid=1;
Calvalid=1;
}
 void IAHP::steptime()
{
assert(Icalvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=Icalvalid=Calvalid=0;
}
void IAHP::getvol()
{
if(!Evalid)
{
extra=new double[2];
extra[1]=tarcal->E();
}
Evalid=1;
}
double IAHP::Icalcium()
{
double v,ca,Ical;
v=tarcal->E();
ca=tarcal->Cain();

Ical=-parameter[0]*(ca/1.+ca)*(v-parameter[1]);
return(Ical);
} 

double IAHP::Calcium()
{
double v,ca,Icala;
v=tarcal->E();
ca=tarcal->Cain();
Icala=-parameter[2]*(1./(1+exp(-(v+25)/2.5)))*(v-parameter[3]);

return(Icala);
}

/*
x[0]= m activation variable
x[1]=h inactivation avriable
parameter[0]= conductance for calcium current in neuron
parameter[1]=Calcium conc conductance term in calcium dynamics
*/
