
#include "modelIzhineuron.h"
#include "math.h"
#include "iostream.h"


inline double MinVal(double a,double b)
{
	double Min;
	if (a<=b)
	Min=a;
	else
	Min=b;
	return(Min);
	
	}

void
modelIzhineuron (double t, double *x, double *dx, double *parameter,
		 double *extra)
{
  double  a, b, c, d,  vpeak, Idc, k,vr,vt;
  double Cap;
  double U;
  Cap = parameter[0];
  a = parameter[1];
  b = parameter[2];
  c = parameter[3];
  d = parameter[4];
  vpeak = parameter[5];
  Idc = parameter[6];
  k = parameter[7];
  vr = parameter[8];
  vt = parameter[9];
  int type=int(parameter[10]);
  
  if (type==1)
  {
   if (x[0]>=vpeak)
   {
    x[0]=c;
    x[1]=x[1]+d;
    }
  }   
    
    if (type==2)
    {
    	if (x[0]>=vpeak-0.1*x[1])
    	{
    		x[0]=-53+.04*x[1];
    		x[1]=MinVal(x[1]+20,670);
    		}
    }
    
    if (type==3)
    {
    	if (x[0]>=vpeak)
    	x[0]=c;
    	if (x[0]<-55)
    	U=0;
    	else
    	U=0.025*(x[0]+55)*(x[0]+55)*(x[0]+55);
    	}
    
  //cout<<x[0]<<endl;  
  dx[0] = (k * (x[0] - vr) * (x[0] - vt) - x[1] + parameter[6] + extra[1])/Cap;
  if (type!=3)
{  dx[1] = a * (b * (x[0] - vr) - x[1]);}
  if (type==3)
{
  dx[1]=0.2*(U-x[1]);
} 
  
}

