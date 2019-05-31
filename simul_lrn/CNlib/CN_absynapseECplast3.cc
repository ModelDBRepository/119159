/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_ABSYNAPSEECPLAST3_CC
#define CN_ABSYNAPSEECPLAST3_CC

#include "CN_absynapse.cc"

#define TINFIN 1e20

absynapseECplast3::absynapseECplast3(neuron *insource, neuron *intarget,
				   double inksyn, double inEsyn, double inEpre,
				   double inasyn, double inbsyn,
				   double inVslope, 
				   double inlrnampl, double indelayT):
  absynapse(insource, intarget, inksyn, inEsyn, inEpre, inasyn,
	    inbsyn, inVslope,
	    ABECPLAST3IVARNO, ABECPLAST3PNO, ABECPLAST3)
{
  p[6]= inlrnampl;
  p[7]= indelayT;
  nextT= TINFIN;
  synapse_change= 0;
} 

absynapseECplast3::~absynapseECplast3()
{
}

void absynapseECplast3::update_gsyn(double *x)
{
  static double chng; 
  if (x[0] > nextT) {  // x[0] is the time ...
    p[0]+= nextdg;
    if (p[0] < 0.0) p[0]= 0.0;
    if (!chngTq.empty()) {
      nextT= chngTq.pop();
      nextdg= dgq.pop();
    }
    else {
      nextT= TINFIN;
    }
    synapse_change= 1;
  }
  else synapse_change= 0;

  if ((source->start_spiking) || (target->start_spiking)) {
    if ((source->spike_time > 0.0) && (target->spike_time > 0.0)) {
      double tau= target->spike_time-source->spike_time;
      chng= STDP_func(tau);
      if (chng != 0.0) {
	if (nextT == TINFIN) {
	  nextT= tG+p[7];
	  nextdg= chng;
	}
	else {
	  chngTq.push(tG+p[7]);
	  dgq.push(chng);
	}
      }
    }
  }
}

double absynapseECplast3::STDP_func(double t)
{
  if (t > 0.0) {  
    return p[6]*pow(t, 10.0)*exp(-abs(t))*5.14e-9;
    // amplitude is fit to data ...
  }
  else {
    return -p[6]*pow(t, 10.0)*exp(-abs(t))*5.14e-9;
  }
}

// end of class implementation

#undef TINFIN

#endif


