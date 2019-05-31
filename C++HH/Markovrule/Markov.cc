#include <iostream.h>
#include<assert.h>
#include "modelMarkov.h"
#include"../rk4.h"
#include "Markov.h"
Markov::Markov(){
gvalid=0;
int i,Dimno,Parno;
Dimno=Markov_DIM;
Parno=Markov_PAR;
xnew=new double[Dimno];
extra=new double[2];
x=new double[Parno];
dx=new double[Dimno];

parameter=new double[Parno];
x[0]=.01;
x[1]=.01;
x[2]=0.75;
x[3]=.25;
x[4]=0.0;
parameter[0]=10.5;//L term
parameter[1]=4.75;//M term
parameter[2]=1; //a term
parameter[3]=1.5; //b term
parameter[4]=6.7;//epsilon_p term
parameter[5]=13.5;//epsilon_d term
parameter[6]=4; //n term
parameter[7]=8.; //taup term
parameter[8]=30.;//taud term
parameter[9]=0.0;//base calcium level
parameter[11]=1;//alpha term characterising high and high lock distribution



}
Markov::Markov(double *a,double *b)
{
gvalid==0;
int i,Dimno,Parno;
Dimno=Markov_DIM;
Parno=Markov_PAR;
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
Markov::~Markov()
{
delete x;
delete dx;
delete parameter;
delete extra;
delete xnew;
}
void Markov::update(double timestep)
{
int Dimno;
Dimno=Markov_DIM;
step(timestep,x,dx,parameter,extra,xnew,Dimno,&t,&modelMarkov);
gvalid=1;
}
 void Markov::steptime()
{
int Dimno;
Dimno=Markov_DIM; 
assert(gvalid);
for(int i=0;i<Dimno;i++)
{
x[i]=xnew[i];

}
gvalid=0;
}
void Markov::getcal(double &ca)
{
if(!gvalid)
{
extra=new double[2];
extra[1]=ca;

//cout <<extra[1]<<endl;
}
gvalid=1;
}
void Markov::basestrength(double a,double b,double c)
{
parameter[12]=a;
parameter[13]=b;
parameter[14]=c;
}

/*

x[0]=P term
x[1]=D term
x[2]=deltag :effective change in conductance term

*/
