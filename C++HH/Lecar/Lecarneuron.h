#ifndef _Lecarneuron_H_
#define _Lecarneuron_H_

#define QN_DIM 5
#define QN_PAR 15

class Lecarneuron:public Neuron{

public:
double t,*x,*dx,*parameter,*extra,*xnew,spike;
Lecarneuron():Neuron(QN_DIM,QN_PAR){}
Lecarneuron(int ,int):Neuron(QN_DIM,QN_PAR)
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
Lecarneuron(double, double); //Constructor of Lecarneuron with Default Parameter values
Lecarneuron(double *a,double *b); //HH neuron with parameter values specified
~Lecarneuron()
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
double sp(){}
double u(){return x[1];}
void operator = (Lecarneuron &test);
friend ostream & operator <<(ostream &, Lecarneuron &);
};

#endif
