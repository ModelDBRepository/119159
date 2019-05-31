#include "iostream.h"
#include "../Neuron.h"
#include "../Synapse.h"
#include"../rk4.h"
#include "Insynapse.h"

Insynapse::Insynapse():Synapse(){}
Insynapse::Insynapse(Neuron *obj2,double a):Synapse((Neuron *) NULL, obj2, I_DIM, I_PAR)
{

x=new double[2];
x[1]=a;
}
Insynapse::~Insynapse()
{
delete []x;
}
void Insynapse::setvalues(Neuron *obj2,double a)
{
source=(Neuron *)NULL;
target=obj2;
target->den.add(this);
x=new double[2];
x[1]=a;
}

void Insynapse::set_I(double Idc)
{
x[1]=Idc;
}
double Insynapse::Isyn()
{
double Isync;
Isync=x[1];
return(Isync);
}

void Insynapse::operator =(Insynapse &b)
{
x=new double[2];
for (int i=0;i<2;i++)
x[i]=b.x[i];
}
/*
void Insynapse::operator <<(ostream & dout, Insynapse &b)
{
dout<<b.x[1];
}
*/
