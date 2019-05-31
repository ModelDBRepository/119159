#ifndef _Rose_H_
#define _Rose_H_

#define QN_DIM 4
#define QN_PAR 8

class Rose:public Neuron{

public:
double t,*x,*dx,*parameter,*extra,*xnew,spike;
Rose():Neuron(QN_DIM,QN_PAR){}
Rose(int ,int):Neuron(QN_DIM,QN_PAR)
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
Rose(double, double); //Constructor of Rose with Default Parameter values
Rose(double *a,double *b); //HH neuron with parameter values specified
~Rose()
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
void operator = (Rose &test);
friend ostream & operator <<(ostream &, Rose &);
};

#endif
