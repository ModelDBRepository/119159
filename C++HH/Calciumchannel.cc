#include "Neuron.h"
#include "Calciumchannel.h"
Calciumchannel::Calciumchannel(Neuron * insource,Neuron * intarget,int dno,int pno)
{
Evalid=0;Icalvalid=0,Calvalid=0;;
sourcecal=insource;
if(sourcecal!=NULL) sourcecal->calsource.add(this);
tarcal=intarget;
if(tarcal!=NULL) tarcal->cal.add(this);

Dimno=dno;
Parno=pno;
}
double Calciumchannel::Icalcium()
{return 0;}
double Calciumchannel::Calcium()
{return 0;}
Calciumchannel::~Calciumchannel()
{
if (sourcecal != NULL) sourcecal->calsource.del(this);
  if (tarcal != NULL) tarcal->cal.del(this);

}
