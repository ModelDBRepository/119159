#include <iostream.h>
#include<assert.h>
#include <math.h>
#include "modelPD.h"
#include"../rk4.h"
#include "PD.h"
PD::PD(){
gvalid=0;
int Dimno,Parno;
Dimno=PD_DIM;
Parno=PD_PAR;
xnew=new double[Dimno];
extra=new double[6];
x=new double[Parno];
dx=new double[Dimno];

for (int i=0;i<6;i++)
extra[i]=0;
parameter=new double[Parno];
x[0]=.1;
x[1]=.15;
x[2]=0.0;
parameter[0]=4;//n term
parameter[1]=1;//ap term
parameter[2]=.098;//bp term
parameter[3]=1.25; //ad term
parameter[4]=.035; //bd term
parameter[5]=1.;//a term
parameter[6]=1;//b
}
PD::PD(double *a,double *b)
{
gvalid=0;
int i,Dimno,Parno;
Dimno=PD_DIM;
Parno=PD_PAR;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];
parameter=new double[Parno];

for(i=0;i<Dimno;i++)
x[i]=a[i];
for(i=0;i<Parno;i++) 
parameter[i]=b[i];  
}                     
PD::~PD()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}
void PD::update(double timestep)
{
int Dimno;
Dimno=PD_DIM;
step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelPD);
gvalid=1;
}
 void PD::steptime()
{
int Dimno;
Dimno=PD_DIM; 
assert(gvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];
}
gvalid=0;
}
void PD::getcal(double &vpre,double &vpost)
{
if(!gvalid)
{
extra=new double[5];
extra[1]=vpre;
extra[2]=vpost;
//cout <<extra[1]<<endl;
}
gvalid=1;
}

void PD::getdiff(double &diff)
{
if (fabs(diff)>0 && fabs(diff)<3.5)
extra[3]=0;
else
extra[3]=1;
}

void PD::operator =(PD &b)
{
t=b.t;

for (int i=0;i<PD_DIM;i++)
{
x[i]=b.x[i];
dx[i]=b.dx[i];
xnew[i]=b.xnew[i];
}
for (int i=0;i<6;i++)
extra[i]=b.extra[i];

for(int i=0;i<PD_PAR;i++)
parameter[i]=b.parameter[i];
}


/*

x[0]=P term
x[1]=D term
x[2]=deltag :effective change in conductance term

*/
