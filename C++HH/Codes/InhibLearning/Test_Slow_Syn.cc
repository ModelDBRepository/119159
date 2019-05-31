#include "define.h"#include "Neur_Class.h"#include "Corr_Stat.h"#include <time.h>intmain(int ac, char **av){	time_t          start, end1, end2;	double          dG, dT, gmutualoldpre, gmutualnewpre, gmutualoldpost,	                gmutualnewpost;	int             maxsz, minsz;	int              maxold=0, maxnew, minold=0, minnew;	time(&start);	int32           seed = time(NULL);	TRandomMersenne rg(seed);	string          s1, s2, s3, s4, s5;	if (ac <= 14) {		puts("NN");		puts("endtime");		puts("Iinput external current I_ext Pre neuron");		puts("Baseline current which is kept fixed, to decide the regime Idc post neuron");		puts("self inhibition strength");		puts("mutual inhibition strength");		puts("decay constant for inhibition");		puts("do you want the time series data");		puts("noise level");		puts("1 if Learning present, O otherwise");		puts("output file name");		puts("Self inhibition strength for pre neuron: Heterogeneity through gs");		puts("Tau Rise time");
		puts("Spike time output file name");
		puts("Iteration Number for Noise Study");		exit(0);	}	int             NN = atoi(av[1]);	double          endtime = atof(av[2]);	double          I_ext = atof(av[3]);	double          Idc = atof(av[4]);//	int             Synapse_Type = atoi(av[5]);	double          gself = atof(av[5]);	double          gmutual = atof(av[6]);	double          taus = atof(av[7]);	int             scale = atoi(av[8]);	double          noise = atof(av[9]);	int             Learnpresent = atoi(av[10]);	string          out = av[11];	double          gselfnew = atof(av[12]);	double          tauR = atof(av[13]);
    int NoiseNum=atoi(av[14]);	double          tauD = taus;	double          tim=0.0, timestep = epsilon, Threshold = 0;	double gelec=0.0;	vector < double >x1, x2, w, *T;	vector < double >x1p, x2p, wp, *TP;	T = new vector < double >[3000];	TP = new vector < double >[3000];		int            *bol, *bolp;	bol = new int   [NN + 10];	bolp = new int  [NN + 10];		for (int i = 0; i < NN; i++) {		x1.push_back(0.);		x2.push_back(0.);		w.push_back(0.);		x1p.push_back(0.);		x2p.push_back(0.);		wp.push_back(0.);	}	ofstream        Outfile, Outfile1, Outfile2;Outfile.open(out.c_str(), ios::app);	s1 = "Timeseries";	s2 = DoubleToStdStr(I_ext);	s3 = DoubleToStdStr(noise);	s4 = IntToString(Learnpresent);	s1 = s1 + "_" + s2 + "_" + s3 + "_" + s4 + ".dat";	s5 = "Rates";	s5 = s5 + "_" + s2 + "_" + s3 + "_" + s4 + ".dat";	s5 = "RatesNetwork/" + s5;//Outfile2.open((char *) s5.c_str(), ios: :out);	HHneuron      **neurpre, **neurpost;	neurpre = new HHneuron *[NN];	neurpost = new HHneuron *[NN];	for (int i = 0; i < NN; i++) {		neurpre[i] = new HHneuron(.1, 1);		neurpre[i]->parameter[11] = -65;		neurpre[i]->parameter[7] = 0.;		neurpre[i]->x[0] = -65;		neurpost[i] = new HHneuron(.1, 1);		neurpost[i]->parameter[11] = -65;		neurpost[i]->parameter[7] = 0;		neurpost[i]->x[0] = -55.23;	}	Insynapse     **Inputpre, **Inputpost;	Inputpre = new Insynapse *[NN];	Inputpost = new Insynapse *[NN];	for (int i = 0; i < NN; i++) {		Inputpre[i] = new Insynapse(neurpre[i], 0.0);		Inputpost[i] = new Insynapse(neurpost[i], 0.0);	}	Couplingcurrent ***coupAB, ***coupBA;	coupAB = new Couplingcurrent **[NN];	coupBA = new Couplingcurrent **[NN];	for (int i = 0; i < NN; i++) {		coupAB[i] = new Couplingcurrent *[NN];		coupBA[i] = new Couplingcurrent *[NN];		for (int j = 0; j < NN; j++) {			coupAB[i][j] =				new Couplingcurrent(neurpre[i], neurpost[j], gelec / NN);			coupBA[i][j] =				new Couplingcurrent(neurpost[i], neurpre[j], gelec / NN);		}	}	TwoDsynapse  ***slowself, ***slowselfA, ***slowmutual, ***slowmutualA;	slowself = new TwoDsynapse **[NN];	slowselfA = new TwoDsynapse **[NN];	slowmutual = new TwoDsynapse **[NN];	slowmutualA = new TwoDsynapse **[NN];	for (int i = 0; i < NN; i++) {		slowself[i] = new TwoDsynapse *[NN];		slowselfA[i] = new TwoDsynapse *[NN];		slowmutual[i] = new TwoDsynapse *[NN];		slowmutualA[i] = new TwoDsynapse *[NN];		slowself[i][i] = new TwoDsynapse(neurpre[i], neurpre[i], gselfnew, 0., 0);		slowselfA[i][i] = new TwoDsynapse(neurpost[i], neurpost[i], gself, 0., 0);		slowself[i][i]->parameter[1] = -82;		slowselfA[i][i]->parameter[1] = -82;		slowself[i][i]->parameter[2] = tauR;		slowself[i][i]->parameter[3] = tauD;		slowselfA[i][i]->parameter[2] = tauR;		slowselfA[i][i]->parameter[3] = tauD;		for (int j = 0; j < NN; j++) {			slowmutual[i][j] = new TwoDsynapse(neurpre[i], neurpost[j], 1.0*gmutual / NN, 0., 0);			slowmutualA[i][j] = new TwoDsynapse(neurpost[i], neurpre[j],  gmutual / NN, 0., 0);			slowmutual[i][j]->parameter[1] = -82;			slowmutualA[i][j]->parameter[1] = -82;			slowmutual[i][j]->parameter[2] = tauR;			slowmutual[i][j]->parameter[3] = tauD;			slowmutualA[i][j]->parameter[2] = tauR;			slowmutualA[i][j]->parameter[3] = tauD;		}	}	int             count = 0;	double         *curr_pre, *curr_post, sumprecur=0, sumpostcur=0, Het;	curr_pre = new double[NN + 10];	curr_post = new double[NN + 10];	for (int i = 0; i < NN; i++) {		curr_pre[i] = I_ext + 0.0 * rg.Random();		curr_post[i] = Idc + 0.0 * rg.Random();		sumprecur += curr_pre[i];		sumpostcur += curr_post[i];	}	Het = 100 * (sumprecur - sumpostcur) / (sumprecur + sumpostcur);	/*******************parameters for the TwoDsynapse to work********************/	double         *ref_time_pre, *ref_time_post;	double         *vnew_pre, *vnew_post, *vold_pre, *vold_post;	double         *diff_pre, *diff_post;	int            *bol_pre, *bol_post, *bol_postpre, *bol_prepost;	ref_time_pre = new double[NN];	ref_time_post = new double[NN];	vnew_pre = new double[NN];	vnew_post = new double[NN];	vold_pre = new double[NN];	vold_post = new double[NN];	diff_pre = new double[NN];	diff_post = new double[NN];	bol_pre = new int[NN];	bol_post = new int[NN];	bol_prepost = new int[NN];	bol_postpre = new int[NN];	for (int i = 0; i < NN; i++) {		vnew_pre[i] = 0;		vnew_post[i] = 0;		vold_pre[i] = neurpre[i]->x[0];		vold_post[i] = neurpost[i]->x[0];		bol_pre[i] = bol_post[i] = 0;		bol_prepost[i] = bol_postpre[i] = 0;	}	/*****************************************************************************/	gmutualoldpre = gmutualoldpost = gmutual;

double Hetgself=100*(gselfnew-gself)/(gselfnew+gself);	
	int szpreold=0,szpostold=0;
	double sumcurr=0,sumcurrpost=0;
	vector<double> PreCurr,PostCurr;
	double RuleNoise=0;
	
	while (tim < endtime) {		for (int i = 0; i < NN; i++) {			Inputpre[i]->set_I(0);			Inputpost[i]->set_I(0);		}		if (tim > 199 && tim < endtime) {			for (int i = 0; i < NN; i++) {				Inputpre[i]->set_I(curr_pre[i]);				Inputpost[i]->set_I(curr_post[i]);			}		}		/*if (tim > 199.25 && tim < 200.) {			for (int i = 0; i < NN; i++) {				Inputpre[i]->set_I(50.0);				Inputpost[i]->set_I(50.0);			}		}*/		if (noise != 0.0) {			for (int i = 0; i < NN; i++) {				do {					x1[i] = 2.0 * rg.Random() - 1.0;					x2[i] = 2.0 * rg.Random() - 1.0;					w[i] = x1[i] * x1[i] + x2[i] * x2[i];				}				while (w[i] >= 1.0);				do {					x1p[i] = 2.0 * rg.Random() - 1.0;					x2p[i] = 2.0 * rg.Random() - 1.0;					wp[i] = x1p[i] * x1p[i] + x2p[i] * x2p[i];				}				while (wp[i] >= 1.0);				neurpre[i]->AddNoise(noise, 0., 1., x1[i], x2[i], w[i],						     timestep);				neurpost[i]->AddNoise(noise, 0., 1., x1p[i], x2p[i], wp[i],						      timestep);			}		}				for (int i = 0; i < NN; i++) {			HH_Run(&tim, neurpre[i], timestep);			HH_Run(&tim, neurpost[i], timestep);			TwoD_Run(&tim,slowself[i][i],neurpre[i]->x[0], timestep);			TwoD_Run(&tim, slowselfA[i][i],neurpost[i]->x[0], timestep);			for (int j = 0; j < NN; j++) {				TwoD_Run(&tim,slowmutual[i][j],neurpre[i]->x[0], timestep);				TwoD_Run(&tim, slowmutualA[i][j],neurpost[i]->x[0], timestep);			}		}		if (tim > 200) {			for (int i = 0; i < NN; i++) {				spiketimes(&tim, &neurpre[i]->x[0], Threshold, bol[i], T[i]);				spiketimes(&tim, &neurpost[i]->x[0], Threshold, bolp[i],					   TP[i]);			}		}
		
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

				
				/* Learning Rule On Line */		if (Learnpresent == 1 ) {			dG = 0;
			
						for (int u = 0; u < NN; u++)				for (int v = 0; v < NN; v++) {			
			
					maxsz = max(T[u].size(), TP[v].size());					minsz = min(T[u].size(), TP[v].size());					maxnew = maxsz;					minnew = minsz;					if (maxnew != maxold || minnew != minold) {                          // if(TP[v].size()>1)            			
            			if (TP[v].size()>0 && T[u].size()>0)				
                         dT = TP[v][TP[v].size() - 1] - T[u][T[u].size() - 1];
                         else
                         dT=-1000;						//RuleNoise=rg.Random();
				//	dG = Rule(dT, 10.0, 1.0, 10.0, 1.0, .01);
						dG=Rule(dT,5,.25,5,.25,.015);
											}					gmutualnewpre = gmutualoldpre + dG;					if (gmutualnewpre < 0)						gmutualnewpre = 0;					gmutualnewpost = gmutualoldpost - dG;					if (gmutualnewpost < 0)						gmutualnewpost = 0;					slowmutual[u][v]->parameter[0] = 1.0*gmutualnewpre;					slowmutualA[u][v]->parameter[0] = gmutualnewpost;
								}			maxold = maxnew;			minold = minnew;		} else;		gmutualoldpre = 1.0*gmutualnewpre;		gmutualoldpost = gmutualnewpost;		count += 1;
		
		if (scale==2)
		{
			if (tim>0)
	Outfile<<tim<<" "<<neurpre[0]->x[0]<<" "<<neurpost[0]->x[0]<<" "<<gmutualnewpost<<" "<<gmutualnewpre<<endl;		}
				tim += epsilon;	} 

 
cout<<tim<<" "<<neurpre[0]->x[0]<<" "<<neurpost[0]->x[0]<<" "<<gmutualnewpre<<" "<<gmutualnewpost<<" "<<slowmutualA[0][0]->x[0]<<endl;
         vector < double >Periodpre, Periodpost;	double          Avgperiodpre, Avgperiodpost;	double          sumpree = 0, sumposte = 0;	double        **kura, sumkura = 0, avgkura;	kura = new double *[NN + 10];	for (int j = 0; j < NN + 10; j++)		kura[j] = new double[NN + 10];	for (int i = 0; i < NN; i++) {		Periodpre.push_back(ISI(T[i]));		Periodpost.push_back(ISI(TP[i]));	}
	for (int i = 0; i < NN; i++) {		sumpree += Periodpre[i];		sumposte += Periodpost[i];	}	Avgperiodpre = sumpree / NN;	Avgperiodpost = sumposte / NN;
	
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
Outfile1.close();*/	    	for (int i = 0; i < NN; i++) {		int             findpre = isnan(Periodpre[i]);		int             findpreinf = isinf(1. / Periodpre[i]);		for (int j = 0; j < NN; j++) {			int             findpost = isnan(Periodpost[j]);			if (findpre != 1 && findpost != 1 && findpreinf != 1			    && Periodpost[j] != 0.) {				kura[i][j] = Kuramoto(T[i], TP[j], Avgperiodpost);				sumkura += kura[i][j];			}		}	}	avgkura = sumkura / (NN * NN);	int             MaxTime, MinTime;	MaxTime = max(T[0].size(), TP[0].size());	MinTime = min(T[0].size(), TP[0].size());	//cout << MaxTime << " " << MinTime << endl;	double          diffspike = MaxTime - 2 * MinTime;
	//Outfile<<I_ext<<" "<<Avgperiodpre<<endl;	
	if (diffspike < 100 && MinTime > 100) {		cout<< NoiseNum<<" "<<gmutual<<" "<<I_ext << " " << Het << " " <<gmutualnewpre << " "<<gmutualnewpost<<" " <<avgkura<<" "<< Avgperiodpre << " " <<Avgperiodpost << " "  << Avgperiodpost / Avgperiodpre << endl;	} 			double         *rateA, *rateB;	rateA = new double[NN + 10];	rateB = new double[NN + 10];	Outfile.close();	Outfile1.close();	time(&end1);double diff = difftime (end1, start);cout<<"Computation time for "<< I_ext<<" is "<<diff<<" sec"<<endl;	delete         *neurpre; delete *neurpost;delete *Inputpre;delete *Inputpost;delete **coupAB;	   delete           **coupBA; delete **slowself; delete **slowselfA; delete **slowmutual;	      delete        **slowmutualA;	delete[] neurpre;	delete[] neurpost;	delete[] Inputpre;	delete[] Inputpost;	delete[] coupAB;	delete[] coupBA;	delete[] slowself;	delete[] slowselfA;	delete[] slowmutual;	delete[] slowmutualA;	delete[] ref_time_pre;	delete[] ref_time_post;	delete[] vnew_pre;	delete[] vnew_post;	delete[] vold_pre;	delete[] vold_post;	delete[] diff_pre;	delete[] diff_post;	delete[] bol_pre;	delete[] bol_post;	delete[] bol_postpre;	delete[] bol_prepost;        delete[] T;        delete []TP;        delete []bol;        delete []bolp;	return 0;}