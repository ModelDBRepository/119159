#include "iostream.h"
#include "modelMinneuron.h"
#include "../rk4.h"
#include "../Neuron.h"
#include "Minneuron.h"

Minneuron::Minneuron(double initial,double iapp):Neuron(QN_DIM,QN_PAR)
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
parameter[0]=20.0;  //gna
parameter[1]=60;  //vna
parameter[2]=10.0;   //gk
parameter[3]=-90.0; //vk
parameter[4]=8.0;  //gl
parameter[5]=-80.0;//vl
parameter[6]=iapp;
parameter[7]=-20;//vm1
parameter[8]=15;//km1
parameter[9]=-25;//vn1
parameter[10]=5;//kn1
parameter[11]=1;//CM
parameter[12]=1;//tn
}
Minneuron::Minneuron(double *a,double *b):Neuron(QN_DIM,QN_PAR)

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
void Minneuron::setvalues(double *a,double *b)
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


void Minneuron::update(double timestep)
{
step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelMinneuron);
Evalid=1;
}
void Minneuron::steptime()
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
void Minneuron::Isynintegrate()
{
if(!Isynvalid)
{
extra[1]=0;
forall(den_it) extra[1]+= den_it->c_value()->Isyn();
}
Isynvalid=1;
}


void Minneuron::operator =(Minneuron &b)
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
ostream &operator <<(ostream & dout, Minneuron &b) //defining '<<' operator for class
                                               //output 
{
dout<<b.x[0];
return(dout);
}
