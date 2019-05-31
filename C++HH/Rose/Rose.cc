#include "iostream.h"
#include "modelRose.h"
#include "../rk4.h"
#include "../Neuron.h"
#include "Rose.h"

Rose::Rose(double initial,double iapp):Neuron(QN_DIM,QN_PAR)
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
x[1]=2;
x[2]=4;
parameter[0]=3.0;  //A
parameter[1]=5;  //B
parameter[2]=-1.6;   //C
parameter[3]=0.0021; //R
parameter[4]=4.0;  //S
parameter[5]=3.281;//I
}


Rose::Rose(double *a,double *b):Neuron(QN_DIM,QN_PAR)
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

void Rose::setvalues(double *a,double *b)
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


void Rose::update(double timestep)
{
step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelRose);
Evalid=1;
}
void Rose::steptime()
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

void Rose::Isynintegrate()
{
if(!Isynvalid)
{
extra[1]=0;
forall(den_it) extra[1]+= den_it->c_value()->Isyn();
}
Isynvalid=1;
}


void Rose::operator =(Rose &b)
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
ostream &operator <<(ostream & dout, Rose &b) //defining '<<' operator for class
                                               //output 
{
dout<<b.x[0];
return(dout);
}
