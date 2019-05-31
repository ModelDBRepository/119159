#include "Connectivity.h"
Connectivity::Connectivity(){bol2D=0;bol1D=0;}
Connectivity::Connectivity(int a,int b,double p)
{
	bol2D=1;
	Nrow=a;
	Ncolumn=b;
	P=p; //probability of connectivity
	Twodarray=new double*[Nrow];
	for(int i=0;i<Nrow;i++)
		Twodarray[i]=new double[Ncolumn];
	for (int i=0;i<Nrow;i++)
		for (int j=0;j<Ncolumn;j++)
			Twodarray[i][j]=0;
	
	seed=time(NULL);
}

Connectivity::Connectivity(int a,double p)
{
	int32 seed1;
	int gen;
	bol1D=1;
	
	TRandomMersenne rg1 (seed1);
	gen=rg1.IRandom(0,99);
	P=p;
	srow=a;
	Onedarray=new double[a];
	
	for (int i=0;i<srow;i++)
		Onedarray[i]=0;
	
	seed=time(NULL)+gen;
}

Connectivity::~Connectivity()
{
	//cout<<"Entered the destructor file"<<endl;
	if (bol2D==1)
	{
		for (int i=0;i<Nrow;i++)
			delete Twodarray[i];
		delete[]Twodarray;
	}
	//if (bol1D==1)delete Onedarray;
}

void Connectivity::setRandomOned()
{
	double gen;
	TRandomMersenne rg (seed);
	for(int i=0;i<srow;i++)
	{
		gen=rg.Random();
		if(gen>P)
			Onedarray[i]=1;
		else
			Onedarray[i]=0;
	}
}

void Connectivity::setRandomTwod()
{
	double gen;
	int32 seed1=time(NULL);
	TRandomMersenne rg (seed1);
	for(int i=0;i<Nrow;i++)
		for(int j=0;j<Ncolumn;j++)
		{
			gen=rg.Random();
			if(gen>P)
			{
				Twodarray[i][j]=1;
			}
			else
				Twodarray[i][j]=0;
		}
}


void Connectivity::setRandomTwodSym() //Symmetric connectivity matrix for gap junction coupling
{
	double gen;
	int32 seed1=time(NULL);
	TRandomMersenne rg (seed1);
	for(int i=0;i<Nrow;i++)
	{
		for(int j=0;j<Ncolumn;j++)
		{
			gen=rg.Random();
			if(gen>P)
			{
					Twodarray[i][j]=1;
					Twodarray[j][i]=1;
			}
			
		}
    }
 for (int i=0;i<Nrow;i++)
 for (int j=0;j<Ncolumn;j++)
 {
 if (i==j)
 Twodarray[i][j]=0;
 }

}


void Connectivity::setStaticOned(vector<int>&A)
{
	if (A.size()!=(unsigned int)srow)
	{
		cout<<"Error in one d array setStatic OneD"<<endl;
		exit(1);
	}
	for (int i=0;i<srow;i++)
		if (A[i]==1)
			Onedarray[i]=1;
		else
			Onedarray[i]=0;
}

void Connectivity::setStaticTwod(double **A)
{
	for (int i=0;i<Nrow;i++)
		for (int j=0;j<Ncolumn;j++)
		{
		Twodarray[i][j]=A[i][j];
		
		/*	if (A[i][j]==1)
				Twodarray[i][j]=1;
			else
				Twodarray[i][j]=0;
		*/	
		}
			
}


void Connectivity::Print()
{
	for (int i=0;i<Nrow;i++)
	{
		for (int j=0;j<Ncolumn;j++)
			cout<<Twodarray[i][j]<<" ";
		cout<<endl;
	}
	
}

