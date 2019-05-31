#include "iostream.h"
//#include "../define.h"
#include "../Neuron.h"
#include "../Synapse.h"
#include "../HHneuron/HHneuron.h"
#include "../Izhineuron/Izhineuron.h"
#include "modelSTFDsynapse.h"
#include "STFDsynapse.h"
#include"../rk4.h"

STFDsynapse::STFDsynapse():Synapse(){}
STFDsynapse::STFDsynapse(Neuron *obj1,Neuron *obj2,double a,double b,int s):Synapse(obj1,obj2,A_DIM,A_PAR)
{
Evalid=Isynvalid=0;
//int i;
steptransfer=s;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];
extra[0]=double(steptransfer);
extra[1]=0;
x[0]=1.0;
x[1]=0.0;
x[2]=0.00;
x[3]=0.00;


parameter[0]=a;// The ASE term for the net synaptic current 
parameter[1]=b; //The USE term for the fraction of neurotransmitters available
parameter[2]=3.0; //the inactivation decay time
parameter[3]=800; // the recovery time scale
parameter[4]=500; //facilitation time scale
parameter[5]=0; //reversal potential for excitatory or inhibitory currents 
}

STFDsynapse::STFDsynapse(Neuron *obj1,Neuron *obj2,double *a,double *b,int s):Synapse(obj1,obj2,A_DIM,A_PAR)
{
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
parameter[i]=b[i];  // Ereversal for the synapse,parameter 2 is Vth of STFDsynapse
}                   //parameter 3 is synaptic timescale for synapse  

STFDsynapse::~STFDsynapse()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}

void STFDsynapse::setvalues(Neuron *obj1, Neuron *obj2,double *a,double *b,int s)
{
int i;
steptransfer=s;
Evalid=Isynvalid=0;
source=obj1;
target=obj2;
source->axon.add(this);
target->den.add(this);

Dimno=A_DIM;
Parno=A_PAR;
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
parameter[i]=b[i];
}

void STFDsynapse::update(double timestep)
{
step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelSTFDsynapse);
Isynvalid=1;
}
 void STFDsynapse::steptime()
{
assert(Isynvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=Isynvalid=0;
}
void STFDsynapse::getvol()
{
if(!Evalid)
{
extra[1]=source->E();
}
Evalid=1;
}
double STFDsynapse::Isyn()
{
double Isync;
Isync=-parameter[0]*x[1]*(target->E()-parameter[5]);;
return(Isync);
}


double STFDsynapse::Calcium()
{
double Isync;
Isync=0.0;//-parameter[6]*x[0]*(target->E()-parameter[1]);
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

