
#include "Neuron.h"
#include "Synapse.h"

Neuron::Neuron(int Dno,int Pno)
{
Evalid=0;
Isynvalid=0;
Iionvalid=0;
Icalvalid=0;
Calvalid=0;
Dimno=Dno;Parno=Pno;
den_it= den.iterator();
  axon_it= axon.iterator();
calden_it=calden.iterator();
calaxon_it=calaxon.iterator();
ion_it=ion.iterator();
cal_it=cal.iterator();
calsource_it=calsource.iterator();
}

Neuron::~Neuron()
{
den_it->init();
  while (den_it->step())
    den_it->c_value()->target= NULL;
  delete den_it;
  den_it= axon.iterator();
  while (den_it->step())
    den_it->c_value()->source= NULL;
  delete den_it;
ion_it->init();
while (ion_it->step())
    ion_it->c_value()->tar= NULL;
  delete ion_it;
  cal_it->init();
while (cal_it->step())
    cal_it->c_value()->tarcal= NULL;
  delete cal_it;
  calsource_it->init();
  while (calsource_it->step())
    calsource_it->c_value()->sourcecal= NULL;
  delete calsource_it;
calden_it->init();
while(calden_it->step())
calden_it->c_value()->source=NULL;
delete calden_it;
calaxon_it->init();
while(calaxon_it->step())
calaxon_it->c_value()->source=NULL;
delete calaxon_it;
  
}
