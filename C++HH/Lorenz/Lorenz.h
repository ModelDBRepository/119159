#ifndef _Lorenz_H_
#define _Lorenz_H_

#define QL_DIM 3
#define QL_PAR 3

class Lorenz:public Neuron{

public:
double t,*x,*dx,*parameter,*extra,*xnew,spike;
Lorenz():Neuron(QL_DIM,QL_PAR){}
Lorenz(int ,int):Neuron(QL_DIM,QL_PAR)
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
Lorenz(double); //Constructor of Lorenz with Default Parameter values
Lorenz(double *a,double *b); //HH neuron with parameter values specified
~Lorenz()
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
double sp(){return 0;}
double u(){return x[1];}
void operator = (Lorenz &test);
friend ostream & operator <<(ostream &, Lorenz &);
};

#endif
