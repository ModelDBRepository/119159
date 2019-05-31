#ifndef _Neur_Class_h
#define _Neur_Class_h
#include "define.h"
void HH_Run(double *tim,HHneuron *neur,double timestep)
{
          neur->intime (*tim);
	  neur->Isynintegrate ();
	  neur->Iionintegrate ();
	  neur->Icalintegrate ();
    neur->Calintegrate();
	  neur->update (timestep);
	  neur->steptime ();
}
/*
void HH_Quadratic_Run(double *tim,QuadraticNeuron *neur,double timestep)
{
          neur->intime (*tim);
	  neur->Isynintegrate ();
	  neur->update (timestep);
	  neur->steptime ();
}



void HHmin_Run(double *tim,Minneuron *neur,double timestep)
{
          neur->intime (*tim);
	  neur->Isynintegrate ();
	  neur->update (timestep);
	  neur->steptime ();
}
*/
void AMPA_Run(double *tim,AMPAsynapse *Ampa,double timestep)
{          
          Ampa->getvol ();
	  Ampa->update (timestep);
	  Ampa->steptime ();
}
void NMDA_Run(double *tim,NMDAsynapse *Nmda,double timestep)
{
          Nmda->getvol ();
	  Nmda->update (timestep);
	  Nmda->steptime ();

}

void GABAA_Run(double *tim,GABAAsynapse *Gaba,double timestep)
{
Gaba->getvol ();
	  Gaba->update (timestep);
	  Gaba->steptime ();

}

void TwoD_Run(double *tim,TwoDsynapse *TwoD,double V,double timestep)
{          
          TwoD->getvol (V);
	  TwoD->update (timestep);
	  TwoD->steptime ();
}

void Henrysynapse_Run(double *tim,double *spk,Henrysynapse *Ex,double timestep)
{
          Ex->set_spike (*spk);
	  Ex->getvol ();
	  Ex->update (timestep);
	  Ex->steptime ();

}
/*
void IA_Run(double *tim,IA *ia,double timestep)
{
ia->getvol();
ia->update(timestep);
ia->steptime(); 
}
void IM_Run(double *tim,IM *im,double timestep)
{
im->getvol();
im->update(timestep);
im->steptime(); 
}

void IL_Run(double *tim,IL *il,double timestep)
{
il->getvol();
il->update(timestep);
il->steptime(); 
}

void IAHP_Run(double *tim,IAHP *iahp,double timestep)
{
iahp->getvol();
iahp->update(timestep);
iahp->steptime(); 
}


void PD_Run(double *tim,double &pre,double &post,double &diff,PD *Pd,double timestep)
{

Pd->intime(*tim);
Pd->getcal(pre,post);
Pd->getdiff(diff);
Pd->update(timestep);
Pd->steptime();
}
*/
#endif
