#ifndef _Neuron_H_
#define _Neuron_H_
#include <iostream.h>
#include<math.h>
#include "Synapse.h"
#include "Ionchannel.h"
#include "Calciumchannel.h"
#include "../simul_lrn/datalib/bst_set.h"

class Calciumchannel;
class Ionchannel;
class Synapse;
class Neuron;



class Neuron{
public:
int Dimno,Parno,Evalid,Isynvalid,Iionvalid,Icalvalid,Calvalid;
  bst_set<Synapse *> den;       // list of all incoming synapses
  bst_set<Synapse *> axon;      // list of all neurons projecting to
bst_set<Synapse *> calden; 
bst_set<Synapse *> calaxon;


bst_set<Ionchannel *> ion;
bst_set<Calciumchannel *> cal;
bst_set<Calciumchannel *> calsource;

Neuron(int,int);
virtual ~Neuron();
  bst_iterator<Synapse *> *den_it;
  bst_iterator<Synapse *> *axon_it;
  bst_iterator<Synapse *> *calden_it;
  bst_iterator<Synapse *> *calaxon_it;
    bst_iterator<Ionchannel *> *ion_it;
    bst_iterator<Calciumchannel *> *cal_it;
bst_iterator<Calciumchannel *> *calsource_it;
virtual void Isynintegrate() = 0;
/*virtual void Iionintegrate()=0;
virtual void Icalintegrate()=0;
*/
virtual double E()=0;
virtual double sp()=0;
virtual double Cain()=0;
};
#endif
