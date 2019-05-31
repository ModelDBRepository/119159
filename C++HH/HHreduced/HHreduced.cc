#include "iostream.h"
#include "modelHHreduced.h"
#include "../rk4.h"
#include "../Neuron.h"
#include "HHreduced.h"

HHreduced::HHreduced(double Idc,int type):Neuron(HH_DIM,HH_PAR)
{
double ENa,Ek,El,gNa,gK,gL,CM,Co;
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

parameter[12]=type; //type of neuron model is characterised. 
                 //==1 => Traub model,==2 implies Henrymodel ==>3 or others
		 //imply some other model  

if (parameter[12]==1)
{gNa=215;gK=43;gL=.813;ENa=50;Ek=-95;El=-64;CM=1;Co=1;}
else if(parameter[12]==2)
{gNa=20;gK=6.2;gL=.03;ENa=50.;Ek=-77.;El=-49.4;CM=1;Co=1;}
else;
x[0]=-85;
x[1]=.034;
x[2]=.33;
x[3]=.03;
x[4]=0;

parameter[0]=CM; //input 0
parameter[1]=ENa;  //.
parameter[2]=gNa;  //.
parameter[3]=Ek;   //.
parameter[4]=gK;   //.
parameter[5]=El;   //.
parameter[6]=gL; //input 6
parameter[7]=Idc; //input 7 is the fixed Idc tothe neuron
parameter[8]=Co; //intracellular equilibrium conc of Calcium
parameter[9]=30; //Calcium decay constant
parameter[10]=1;//scaling for Calcium conc through various Calcium sources
parameter[11]=-48; //The Vth term in Traub model used to adjust the threshold of spiking

}
HHreduced::HHreduced(double *a,double *b):Neuron(HH_DIM,HH_PAR)

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
void HHreduced::setvalues(double *a,double *b)
{
Evalid=Iionvalid=Isynvalid=0;
int i;
Dimno=HH_DIM;
Parno=HH_PAR;
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


void HHreduced::update(double timestep)
{
step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelHHreduced);
Evalid=1;
}
void HHreduced::steptime()
{
assert(Evalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
Evalid=0;
Isynvalid=0;
Iionvalid=0;
Icalvalid=0;
Calvalid=0;
}
void HHreduced::Isynintegrate()
{
if(!Isynvalid)
{
extra[1]=0;
forall(den_it) extra[1]+= den_it->c_value()->Isyn();
}
Isynvalid=1;
}
void HHreduced::Iionintegrate()
{
if(!Iionvalid)
{
extra[2]=0;
forall(ion_it) extra[2]+= ion_it->c_value()->Iion();
}
Iionvalid=1;
}
void HHreduced::Icalintegrate()
{
if(!Icalvalid)
{
extra[3]=0;
forall(cal_it) extra[3]+= cal_it->c_value()->Icalcium();
}
Icalvalid=1;
}

void HHreduced::Calintegrate()
{
if(!Calvalid)
{
extra[4]=0;
forall(cal_it) extra[4]+= cal_it->c_value()->Calcium();
forall(calden_it)extra[4]+=calden_it->c_value()->Calcium();
//cout<<extra[4]<<endl;
}
Calvalid=1;
}

void HHreduced::operator =(HHreduced &b)
{
t=b.t;
x=new double[HH_DIM];
dx=new double[HH_DIM];
xnew=new double[HH_DIM];
extra=new double[5];
parameter=new double[HH_PAR];
for (int i=0;i<HH_DIM;i++)
{
x[i]=b.x[i];
dx[i]=b.dx[i];
xnew[i]=b.xnew[i];
}
extra[1]=b.extra[1];
for(int i=0;i<HH_PAR;i++)
parameter[i]=b.parameter[i];
}
ostream &operator <<(ostream & dout, HHreduced &b) //defining '<<' operator for class
                                               //output 
{
dout<<b.x[0];
return(dout);
}
