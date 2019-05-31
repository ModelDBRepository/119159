#include "iostream.h"
#include "modelIzhineuron.h"
#include "../rk4.h"
#include "../Neuron.h"
#include "Izhineuron.h"

Izhineuron::Izhineuron(double Cap,double a,double b,double c,double d, double vpeak,double Idc,double k,double vr,double vt,int type):Neuron(Izhi_DIM,Izhi_PAR)
{
int i;
Evalid=0;
Isynvalid=0;
Iionvalid=0;
t=0;
extra=new double[5];
for(i=0;i<5;i++)
extra[i]=0;
xnew=new double[Parno];
parameter=new double[Parno];
x=new double[Dimno];
dx=new double[Dimno];

parameter[0]=Cap; //Capacitance term
//The following numbering is used
/*1=SN,2=AF,3=DLMP,4=DLMIN,5=lMANP,6=RAP,7=RAIN,*/  

x[0]=vr;
x[1]=0;
x[2]=0;
parameter[1]=a;  //.
parameter[2]=b;  //.
parameter[3]=c;   //.
parameter[4]=d;   //.
parameter[5]=vpeak;   //.
//parameter[6]=alpha; //input 6
parameter[6]=Idc; //input 7 is the fixed Idc tothe neuron
parameter[7]=k;
parameter[8]=vr;
parameter[9]=vt;
parameter[10]=double(type);

}
Izhineuron::Izhineuron(double *a,double *b):Neuron(Izhi_DIM,Izhi_PAR)

{
//cout<<"evoked constructor" <<endl;

Evalid=Iionvalid=Isynvalid=0;
int i;
extra=new double[5];
xnew=new double[Dimno];
parameter=new double[Parno];
x=new double[Dimno];
dx=new double[Parno];
t=0;
for(i=0;i<5;i++)
extra[i]=0;
for(i=0;i<Parno;i++)
{
parameter[i]=b[i];
}
for(i=0;i<Dimno;i++)
{
x[i]=a[i];
}
}
void Izhineuron::setvalues(double *a,double *b)
{
Evalid=Iionvalid=Isynvalid=0;
int i;
Dimno=Izhi_DIM;
Parno=Izhi_PAR;
extra=new double[5];
xnew=new double[Dimno];
parameter=new double[Parno];
x=new double[Dimno];
dx=new double[Parno];
t=0;
for(i=0;i<5;i++)
extra[i]=0;
for(i=0;i<Parno;i++)
{
parameter[i]=b[i];
}
for(i=0;i<Dimno;i++)
{
x[i]=a[i];
}

}


void Izhineuron::update(double timestep)
{
step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelIzhineuron);
Evalid=1;
}
void Izhineuron::steptime()
{
assert(Evalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=0;
Isynvalid=0;
Iionvalid=0;

}
void Izhineuron::Isynintegrate()
{
	int count=0;
if(!Isynvalid)
{
extra[1]=0;

forall(den_it){ extra[1]+= den_it->c_value()->Isyn();}
}
Isynvalid=1;
}



void Izhineuron::operator =(Izhineuron &b)
{
t=b.t;

for (int i=0;i<Izhi_DIM;i++)
{
x[i]=b.x[i];
dx[i]=b.dx[i];
xnew[i]=b.xnew[i];
}
extra[1]=b.extra[1];
for(int i=0;i<Izhi_PAR;i++)
parameter[i]=b.parameter[i];
}
ostream &operator <<(ostream & dout, Izhineuron &b) //defining '<<' operator for class
                                               //output 
{
dout<<b.x[0];
return(dout);
}
