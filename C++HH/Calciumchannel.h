#ifndef _Calciumchannel_H_
#define _Calciumchannel_H_

#include "Neuron.h"

class Neuron;
class Calciumchannel;
class Calciumchannel{
public:
int Dimno,Evalid,Icalvalid,Calvalid;
int Parno;
Neuron *sourcecal;
Neuron *tarcal;
Calciumchannel(Neuron *,Neuron *,int ,int);
Calciumchannel(){}
virtual ~Calciumchannel();
virtual double Icalcium();
virtual double Calcium();
};

#endif
