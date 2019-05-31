#include "iostream.h"
#include "../Neuron.h"
#include "../Synapse.h"
#include "../HHneuron/HHneuron.h"
#include"../rk4.h"
#include "Couplingcurrent.h"

Couplingcurrent::Couplingcurrent():Synapse(){}
Couplingcurrent::Couplingcurrent(Neuron *obj1,Neuron *obj2,double a):Synapse(obj1,obj2,I_DIM,I_PAR)
{
source=obj1;
target=obj2;
source->axon.add(this);
target->den.add(this);
parameter=new double[2];
parameter[0]=a;
}
Couplingcurrent::~Couplingcurrent()
{
delete []parameter;
}
void Couplingcurrent::setvalues(Neuron *obj1,Neuron *obj2,double a)
{
source=obj1;
target=obj2;
source->axon.add(this);
target->den.add(this);
parameter=new double[2];
parameter[0]=a;
}


double Couplingcurrent::Isyn()
{
double Isync;
Isync=parameter[0]*(source->E()-target->E());
return(Isync);
}
