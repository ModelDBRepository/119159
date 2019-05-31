#include "iostream.h"
//#include "../define.h"
#include "../Neuron.h"
#include "../Synapse.h"
#include "../HHneuron/HHneuron.h"
#include "../Izhineuron/Izhineuron.h"
#include "modelTwoDsynapse.h"
#include"../rk4.h"
#include "TwoDsynapse.h"
TwoDsynapse::TwoDsynapse():Synapse(){}
TwoDsynapse::TwoDsynapse(Neuron *obj1,Neuron *obj2,double a,double b,int s ):Synapse(obj1,obj2,Two_DIM,Two_PAR)
{
	Dimno=Two_DIM;
Parno=Two_PAR;

Evalid=Isynvalid=0;

steptransfer=s;
xnew=new double[Dimno];
extra=new double[3]; //This was the culprit!!!!!!!
x=new double[Parno]; //I had the problem of Mac compiler and the linux compiler giving out different results
dx=new double[Dimno]; //It turned out I was using extra[2] variable when I had defined only extra=new double[2]
parameter=new double[Parno]; //I had not assigned memory to use exta[2] in the first place!!!!!talathi
extra[0]=double(steptransfer);
extra[1]=0;
extra[2]=0;
x[0]=.00;
parameter[0]=a;
parameter[1]=0;
parameter[2]=.1;//rise time for synapse
parameter[3]=1.5;//fall time for synapse
parameter[4]=2;
parameter[5]=5;
parameter[6]=b;
}

TwoDsynapse::TwoDsynapse(Neuron *obj1,Neuron *obj2,double *a,double *b,int s):Synapse(obj1,obj2,Two_DIM,Two_PAR)
{
	Dimno=Two_DIM;
Parno=Two_PAR;

Evalid=Isynvalid=0;
int i;
steptransfer=s;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
extra[0]=double(steptransfer);
extra[1]=0;
for(i=0;i<Dimno;i++)
x[i]=a[i];
for(i=0;i<Parno;i++) //parameter 0 is the synaptic strength and parameter 1 is
parameter[i]=b[i];  // Ereversal for the synapse,parameter 2 is Vth of TwoDsynapse
}                   //parameter 3 is synaptic timescale for synapse  
TwoDsynapse::~TwoDsynapse()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}
void TwoDsynapse::setvalues(Neuron *obj1, Neuron *obj2,double *a,double *b,int s)
{
int i;
steptransfer=s;
Evalid=Isynvalid=0;
source=obj1;
target=obj2;
source->axon.add(this);
target->den.add(this);

Dimno=Two_DIM;
Parno=Two_PAR;
xnew=new double[Dimno];
extra=new double[3];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
extra[0]=double(steptransfer);
extra[1]=0;
extra[2]=0;
for(i=0;i<Dimno;i++)
x[i]=a[i];
for(i=0;i<Parno;i++) //parameter 0 is the synaptic strength and parameter 1 is
parameter[i]=b[i];
}
void TwoDsynapse::update(double timestep)
{
	Dimno=Two_DIM;
Parno=Two_PAR;

step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelTwoDsynapse);
Isynvalid=1;
}
 void TwoDsynapse::steptime()
{
assert(Isynvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=Isynvalid=0;
}
void TwoDsynapse::getvol(double V)
{
if(!Evalid)
{
//extra[1]=source->E();
extra[1]=V;
}
Evalid=1;
}
double TwoDsynapse::Isyn()
{
double Isync;
Isync=-parameter[0]*x[0]*(target->E()-parameter[1]);
return(Isync);
}


double TwoDsynapse::Calcium()
{
double Isync;
Isync=-parameter[6]*(1.0*x[0])*(target->E()-parameter[1]);
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

*/

