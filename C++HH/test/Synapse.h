#ifndef _Synapse_H_
#define _Synapse_H_

#include "Neuron.h"

class Neuron;
class Synapse;
class Synapse{
public:
int Dimno,Evalid,Isynvalid;
int Parno;

Neuron *source;
Neuron *target;
Synapse(Neuron *,Neuron *,int ,int);
Synapse(){}
virtual ~Synapse();
virtual double Isyn();
virtual double Calcium();
};

#endif
