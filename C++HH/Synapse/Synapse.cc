
#include "Neuron.h"
#include "Synapse.h"
Synapse::Synapse(Neuron *insource,Neuron * intarget,int dno,int pno)
{
Evalid=Isynvalid=0;
source=insource;
if(source!=NULL) {source->axon.add(this);source->calaxon.add(this);}
target=intarget;
if(target!=NULL){ target->den.add(this);target->calden.add(this);}
Dimno=dno;
Parno=pno;
}

Synapse::~Synapse()
{
if (source != NULL) {source->axon.del(this);source->calaxon.del(this);}
  if (target != NULL) {target->den.del(this);target->calden.del(this);}
}
double Synapse::Isyn()
{return 0;}
double Synapse::Calcium()
{return 0;}
