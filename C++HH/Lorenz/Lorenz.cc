#include "iostream.h"
#include "modelLorenz.h"
#include "../rk4.h"
#include "../Neuron.h"
#include "Lorenz.h"

Lorenz::Lorenz(double initial):Neuron(QL_DIM,QL_PAR)
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
x[1]=0.5;
x[2]=0.5;
parameter[0]=28.0;  //rho
parameter[1]=8.0/3.0;  //beta
parameter[2]=10.0;   //sigma

}


Lorenz::Lorenz(double *a,double *b):Neuron(QL_DIM,QL_PAR)
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

void Lorenz::setvalues(double *a,double *b)
{
Evalid=Iionvalid=Isynvalid=0;
int i;
Dimno=QL_DIM;
Parno=QL_PAR;
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


void Lorenz::update(double timestep)
{
step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelLorenz);
Evalid=1;
}
void Lorenz::steptime()
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

void Lorenz::Isynintegrate()
{
if(!Isynvalid)
{
extra[1]=0;
forall(den_it) extra[1]+= den_it->c_value()->Isyn();
}
Isynvalid=1;
}


void Lorenz::operator =(Lorenz &b)
{
t=b.t;

for (int i=0;i<QL_DIM;i++)
{
x[i]=b.x[i];
dx[i]=b.dx[i];
xnew[i]=b.xnew[i];
}
extra[1]=b.extra[1];
for(int i=0;i<QL_PAR;i++)
parameter[i]=b.parameter[i];
}
ostream &operator <<(ostream & dout, Lorenz &b) //defining '<<' operator for class
                                               //output 
{
dout<<b.x[0];
return(dout);
}
