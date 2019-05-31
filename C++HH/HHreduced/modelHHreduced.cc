
#include "modelHHreduced.h"
#include "math.h"
#include "iostream.h"
void modelHHreduced(double t,double *x,double *dx,double *parameter,double *extra)
{
double am,bm,ah,bh,an,bn,v,minf,h;
v=x[0];
if(parameter[12]==1 ) //traub  reduced model type I
{
am=.32*(13-(v-parameter[11]))/(exp((13-(v-parameter[11]))/4.)-1);
bm=.28*((v-parameter[11])-40)/(exp(((v-parameter[11])-40)/5)-1);
an=0.032*(15-(v-parameter[11]))/(exp((15-(v-parameter[11]))/5)-1);
bn=.5*exp((10-(v-parameter[11]))/40);
minf=am/(am+bm);
h=(-x[1]+.8114)/.8154;
}
else if(parameter[12]==2 ) //henry reduced model type II
{
am=-.1*(35+v)/(exp(-(35+v)/10)-1);
bm=4.*exp(-(60+v)/18);
an=-0.01*(50+v)/(exp(-(50+v)/10)-1);
bn=.125*exp(-(60+v)/80);
minf=am/(am+bm);
h=(-x[1]+.7499)/.8314;
}

dx[0]=(parameter[7]+parameter[2]*pow(minf,3)*h*(parameter[1]-x[0])+
parameter[4]*pow(x[1],4)*(parameter[3]-x[0])+parameter[6]*(parameter[5]-x[0])+extra[1]+extra[2]+1.0*extra[3])/parameter[0];

dx[1]=an*(1-x[1])-bn*x[1];
}

/*
Traub model with Vth=-67
gNa=215;gK=43;gL=.813;ENa=50;Ek=-95;El=-64;CM=1;Co=1;

am=-.32*(52+x[0])/(exp(-(x[0]+52)/4)-1);
bm=.28*(25+x[0])/(exp((25+x[0])/5)-1);
ah=.128*exp(-(x[0]+48)/18);
bh=4.0/(1+exp(-(x[0]+25)/5));
an=-.032*(50+x[0])/(1+exp(-(50+x[0])/5));
bn=.5*exp(-(55+x[0])/40);


The reduced model for this system can be obtained as
replace m variable with minf,
and the following relation between h and n variables

n=-.8154h+.8114

Std HH model 
gNa=20;gK=6.2;gL=.03;ENa=50.;Ek=-77.;El=-49.4;CM=1;Co=1;

am=-.1*(35+v)/(exp(-(35+v)/10)-1);
bm=4.*exp(-(60+v)/18);
ah=.07*exp(-(60+v)/20);
bh=1./(exp(-(30+v)/10)+1);
an=-0.01*(50+v)/(exp(-(50+v)/10)-1);
bn=.125*exp(-(60+v)/80);

The reduced model can againbe obtianed by 
replacing m with minf and
the following relation between h and the n variables.

n=-.8314h+.7499





*/

