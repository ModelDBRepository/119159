/*Role of IA currents */

/*
IA currents play role in getting continuous gain function for spiking neuron...
Also in there presence there is significant delay between firing of neuron in
response to a depolarizing dc current pulse
*/

#ifndef _Markovrule_H_
#define _Markov_H_
#define Markov_DIM 6
#define Markov_PAR 15
class Markov{

public:
double t,*x,*dx,*xnew,*parameter,*extra,gvalid;
Markov();
Markov(double *a, double *b);
~Markov();
void update(double timestep);
 void steptime();
void intime(double intim){t=intim;}
void getcal(double &);
double deltag(){return parameter[12]*x[2]+parameter[13]*x[3]+parameter[14]*x[4];}
void basestrength(double, double,double);
};
#endif
