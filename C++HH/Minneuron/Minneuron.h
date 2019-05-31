#ifndef _Minneuron_H_
#define _Minneuron_H_

#define QN_DIM 5
#define QN_PAR 15

class Minneuron:public Neuron{

public:
double t,*x,*dx,*parameter,*extra,*xnew,spike;
Minneuron():Neuron(QN_DIM,QN_PAR){}
Minneuron(int ,int):Neuron(QN_DIM,QN_PAR)
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
Minneuron(double, double); //Constructor of Minneuron with Default Parameter values
Minneuron(double *a,double *b); //HH neuron with parameter values specified
~Minneuron()
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
void operator = (Minneuron &test);
friend ostream & operator <<(ostream &, Minneuron &);
};

#endif
