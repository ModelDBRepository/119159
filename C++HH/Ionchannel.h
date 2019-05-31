#ifndef _Ionchannel_H_
#define _Ionchannel_H_

#include "Neuron.h"

class Neuron;
class Ionchannel;
class Ionchannel{
public:
int Dimno,Evalid,Iionvalid,Icalvalid;
int Parno;
Neuron *tar;
Ionchannel(Neuron *,int ,int);
Ionchannel(){}
virtual ~Ionchannel();
virtual double Iion()=0;
};

#endif
