#ifndef _Izhineuron_H_
#define _Izhineuron_H_

#define Izhi_DIM 5
#define Izhi_PAR 15

class Izhineuron:public Neuron{

public:
double t,*x,*dx,*parameter,*extra,*xnew,spike;
Izhineuron():Neuron(Izhi_DIM,Izhi_PAR){}
Izhineuron(int ,int):Neuron(Izhi_DIM,Izhi_PAR)
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
Izhineuron(double, double,  double, double,double,double,double,double,double,double,int); //Constructor of Izhineuron with Default Parameter values
Izhineuron(double *a,double *b); //HH neuron with parameter values specified
~Izhineuron()
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
double Cain() {return 0.0;}
double sp(){return x[2];}
double u(){return x[1];}
void operator = (Izhineuron &test);
friend ostream & operator <<(ostream &, Izhineuron &);
};

#endif
