
#include "Neuron.h"
#include "Ionchannel.h"
Ionchannel::Ionchannel(Neuron * intarget,int dno,int pno)
{
Evalid=Iionvalid=0;
tar=intarget;
if(tar!=NULL) tar->ion.add(this);


Dimno=dno;
Parno=pno;
}

Ionchannel::~Ionchannel()
{
  if (tar != NULL) tar->ion.del(this);
  
}
