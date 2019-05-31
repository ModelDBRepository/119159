#include "iostream.h"
#include "modelLecarneuron.h"
#include "../rk4.h"
#include "../Neuron.h"
#include "Lecarneuron.h"

Lecarneuron::Lecarneuron(double initial,double iapp):Neuron(QN_DIM,QN_PAR)
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
parameter[0]=4.0;  //gca
parameter[1]=120;  //vca
parameter[2]=8.0;   //gk
parameter[3]=-84.0; //vk
parameter[4]=2.0;  //gl
parameter[5]=-60.0;//vl
parameter[6]=iapp;
parameter[7]=-1.2;//v1
parameter[8]=18;//v2
parameter[9]=12;//v3
parameter[10]=17.4;//v4
parameter[11]=.06667;//phi
parameter[12]=20;//CM

}
Lecarneuron::Lecarneuron(double *a,double *b):Neuron(QN_DIM,QN_PAR)

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
void Lecarneuron::setvalues(double *a,double *b)
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


void Lecarneuron::update(double timestep)
{
step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelLecarneuron);
Evalid=1;
}
void Lecarneuron::steptime()
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
void Lecarneuron::Isynintegrate()
{
if(!Isynvalid)
{
extra[1]=0;
forall(den_it) extra[1]+= den_it->c_value()->Isyn();
}
Isynvalid=1;
}


void Lecarneuron::operator =(Lecarneuron &b)
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
ostream &operator <<(ostream & dout, Lecarneuron &b) //defining '<<' operator for class
                                               //output 
{
dout<<b.x[0];
return(dout);
}
