#include "define.h"
		puts("Spike time output file name");
		puts("Iteration Number for Noise Study");
    int NoiseNum=atoi(av[14]);

double Hetgself=100*(gselfnew-gself)/(gselfnew+gself);
	int szpreold=0,szpostold=0;
	double sumcurr=0,sumcurrpost=0;
	vector<double> PreCurr,PostCurr;
	double RuleNoise=0;
	
	while (tim < endtime) {
		
		/*Getting the total current per cycle entering the neuron*/
		
		if (szpreold!=T[0].size())
		{
			//cout<<"Entered loop "<<T[0].size()-1<<endl;
		    PreCurr.push_back(sumcurr);
		    sumcurr=0;
		}
		sumcurr+=I_ext+slowself[0][0]->Isyn()+slowmutualA[0][0]->Isyn();
		szpreold=T[0].size();
		
				if (szpostold!=TP[0].size())
		{
			//cout<<"Entered loop "<<TP[0].size()-1<<endl;
		    PostCurr.push_back(sumcurrpost);
		    sumcurrpost=0;
		}
		sumcurrpost+=Idc+slowselfA[0][0]->Isyn()+slowmutual[0][0]->Isyn();
		szpostold=TP[0].size();

				
		
			
			
			
					maxsz = max(T[u].size(), TP[v].size());
            			if (TP[v].size()>0 && T[u].size()>0)				
                         dT = TP[v][TP[v].size() - 1] - T[u][T[u].size() - 1];
                         else
                         dT=-1000;
				//	dG = Rule(dT, 10.0, 1.0, 10.0, 1.0, .01);
						dG=Rule(dT,5,.25,5,.25,.015);
						
				
		
		if (scale==2)
		{
			if (tim>0)
	Outfile<<tim<<" "<<neurpre[0]->x[0]<<" "<<neurpost[0]->x[0]<<" "<<gmutualnewpost<<" "<<gmutualnewpre<<endl;
		

 
cout<<tim<<" "<<neurpre[0]->x[0]<<" "<<neurpost[0]->x[0]<<" "<<gmutualnewpre<<" "<<gmutualnewpost<<" "<<slowmutualA[0][0]->x[0]<<endl;
 

	
//cout<<I_ext<<" "<<Avgperiodpre<<" "<<Avgperiodpost<<endl;

/*
Outfile1.open((char *)SpikeTimeFile.c_str(), ios::out);
if (scale!=0)
{
	for (int i=0;i<PreCurr.size();i++)
	Outfile1<<PreCurr[i]<<endl;
	Outfile1<<-1000<<endl;
	for (int i=0;i<PostCurr.size();i++)
	Outfile1<<PostCurr[i]<<endl;
	
for (int i=0;i<T[0].size();i++)
Outfile1<<T[0][i]<<endl;
Outfile1<<-1000<<endl;
for (int i=0;i<TP[0].size();i++)
{Outfile1<<TP[0][i]<<endl;}
}
Outfile1.close();
	//Outfile<<I_ext<<" "<<Avgperiodpre<<endl;
	if (diffspike < 100 && MinTime > 100) {