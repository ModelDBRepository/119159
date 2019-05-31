#include "iostream.h"
#include "modelQuadraticNeuron.h"
#include "../rk4.h"
#include "../Neuron.h"
#include "QuadraticNeuron.h"

QuadraticNeuron::QuadraticNeuron(double initial,double b,double vreset,double vpeak, double k,double c,double vrest):Neuron(QN_DIM,QN_PAR)
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


x[0]=initial;
x[1]=0;
x[2]=0;
parameter[0]=b;  //.
parameter[1]=vreset;  //.
parameter[2]=vpeak;   //.
parameter[3]=k;
parameter[4]=vrest;
parameter[5]=c;
}
QuadraticNeuron::QuadraticNeuron(double *a,double *b):Neuron(QN_DIM,QN_PAR)

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
void QuadraticNeuron::setvalues(double *a,double *b)
{
Evalid=Iionvalid=Isynvalid=0;
int i;
Dimno=QN_DIM;
Parno=QN_PAR;
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


void QuadraticNeuron::update(double timestep)
{
step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelQuadraticNeuron);
Evalid=1;
}
void QuadraticNeuron::steptime()
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
void QuadraticNeuron::Isynintegrate()
{
if(!Isynvalid)
{
extra[1]=0;
forall(den_it) extra[1]+= den_it->c_value()->Isyn();
}
Isynvalid=1;
}


void QuadraticNeuron::operator =(QuadraticNeuron &b)
{
t=b.t;

for (int i=0;i<QN_DIM;i++)
{
x[i]=b.x[i];
dx[i]=b.dx[i];
xnew[i]=b.xnew[i];
}
extra[1]=b.extra[1];
for(int i=0;i<QN_PAR;i++)
parameter[i]=b.parameter[i];
}
ostream &operator <<(ostream & dout, QuadraticNeuron &b) //defining '<<' operator for class
                                               //output 
{
dout<<b.x[0];
return(dout);
}
