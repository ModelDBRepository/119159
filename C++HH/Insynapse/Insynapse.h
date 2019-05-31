#ifndef _Insynapse_H_
#define _Insynapse_H_
#define I_DIM 2
#define I_PAR 2
class Insynapse:public Synapse{

public:
double *x;
Insynapse();
Insynapse(Neuron *,double );
~Insynapse();
void setvalues(Neuron *,double);
void set_I(double);
virtual double Isyn();
void operator =(Insynapse &test);
//friend ostream & operator <<(ostream &,Insynapse &);
};
#endif
