#ifndef _QuadraticNeuron_H_
#define _QuadraticNeuron_H_

#define QN_DIM 5
#define QN_PAR 15

class QuadraticNeuron:public Neuron{

public:
double t,*x,*dx,*parameter,*extra,*xnew,spike;
QuadraticNeuron():Neuron(QN_DIM,QN_PAR){}
QuadraticNeuron(int ,int):Neuron(QN_DIM,QN_PAR)
{
extra=new double[5];
xnew=new double[Dimno];
parameter=new double[Parno];
x=new double[Dimno];
dx=new double[Parno];
t=0;
for(int i=0;i<5;i++)
extra[i]=0;
for(int i=0;i<Parno;i++)
{
parameter[i]=0;
}
for(int i=0;i<Dimno;i++)
{
x[i]=0;
}
}
QuadraticNeuron(double, double,  double, double,double,double,double); //Constructor of QuadraticNeuron with Default Parameter values
QuadraticNeuron(double *a,double *b); //HH neuron with parameter values specified
~QuadraticNeuron()
{
delete []x;delete []dx; delete []parameter,delete []xnew,delete []extra;

}
void setvalues(double *a,double *b);
void update(double timestep);
void steptime();
void intime(double intim){t=intim;}
double time(){return t;}
void Isynintegrate();
double E(){return x[0];}
double Cain(){return 0;}
double sp(){}
double u(){return x[1];}
void operator = (QuadraticNeuron &test);
friend ostream & operator <<(ostream &, QuadraticNeuron &);
};

#endif
