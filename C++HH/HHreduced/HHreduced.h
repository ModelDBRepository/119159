#ifndef _HHreduced_H_
#define _HHreduced_H_

#define HH_DIM 5
#define HH_PAR 13

class HHreduced:public Neuron{

public:
double t,*x,*dx,*parameter,*extra,*xnew;
HHreduced():Neuron(HH_DIM,HH_PAR){}
HHreduced(int ,int):Neuron(HH_DIM,HH_PAR)
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
HHreduced(double Idc,int ); //Constructor of HHreduced with Default Parameter values
HHreduced(double *a,double *b); //HH neuron with parameter values specified
~HHreduced()
{
delete []x;delete []dx; delete []parameter,delete []xnew,delete []extra;

}
void setvalues(double *a,double *b);
void update(double timestep);
void steptime();
void intime(double intim){t=intim;}
double time(){return t;}
void Isynintegrate();
void Iionintegrate();
void Icalintegrate();
void Calintegrate();
double E(){return x[0];}
double sp(){};
double Cain(){return x[4];}
void operator = (HHreduced &test);
friend ostream & operator <<(ostream &, HHreduced &);
};

#endif
