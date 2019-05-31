#ifndef _HHneuron_H_
#define _HHneuron_H_

#define HH_DIM 5
#define HH_PAR 13

class HHneuron:public Neuron{

public:
double t,*x,*dx,*parameter,*extra,*xnew;
HHneuron():Neuron(HH_DIM,HH_PAR){}
HHneuron(int ,int):Neuron(HH_DIM,HH_PAR)
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
HHneuron(double Idc,int ); //Constructor of HHneuron with Default Parameter values
HHneuron(double *a,double *b); //HH neuron with parameter values specified
~HHneuron()
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
double sp(){double test; test=0; return(test);};
double Cain(){return x[4];}
double box_muller(double m, double s, double &x1, double &x2, double &w);
void AddNoise(double scale,double mean,double std,double &x1,double &x2,double &w,double);
void operator = (HHneuron &test);
friend ostream & operator <<(ostream &, HHneuron &);
};

#endif
