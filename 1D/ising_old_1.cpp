#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <vector>
#include <cstdio>
using namespace std;
int lattice[100],trial[100];
double energystrength=1;
int N= 50 ;//Number of lattice sites 
int maxiter=500;// number of runs
double Temperature=5; // Temperature
double B=1/Temperature ; // 1/KbT 
int flag=0;
double siteenergy=1;



double energy( int lattice[])
{   double energy=0,site=0;
    for ( int i=0 ; i<N ; i++ )
    {   if(i==N-1)
          {energy+= double(lattice[i]*lattice[0]);
          site+=double(lattice[i]);}
        else
          {energy+= double(lattice[i]*lattice[i+1]);
          site+=double(lattice[i]);}  
    } 
    site*=(-siteenergy);
    energy*=-energystrength;  
    energy+=site;  
    return energy;
}
void initialize_hot()
{  
    for(int i=0;i<N;i++)
    {    
        if(rand()%10<5)
          {lattice[i]=1;
           trial[i]=1;}
        else
          {lattice[i]=-1;
           trial[i]=-1;}
    }
}

void initialize_cold()
{  
    for(int i=0;i<N;i++)
          {lattice[i]=1;
           trial[i]=1;}

    
}
int main()
{ ofstream data("data.txt");
  initialize_cold();
  double E = 0.0;
  double Etrial=0.0;   
  double minimumenergy=0.0;
  int flag=0;
  for (int t=0;t<maxiter;t++)
    {   
            
        cout<<"\n";
          for(int i=0;i<N;i++)
            { if (lattice[i]==1)
                cout<<"|"<<1<<"|";
              else 
                cout<<"|"<<0<<"|";
            }

      E = energy(lattice);
        // Choosing a random particle
        int choose=rand()%N;
        // flipping the spin of the random particle
        trial[choose]*=-1;
        // Calculating new configuration energy
        Etrial=energy(trial);

        if(Etrial < E) // checking if the trial energy is less than the original energy  
         {for(int i=0;i<N;i++)
              {lattice[i]=trial[i];} }// accepting the new configuration
          

        else { // If the trial energy is greater then 
        double R=double(rand()%10)/10;  // choosing a number between (0,1) 
          double bweight = exp(-(Etrial-E)*B); // calculating the boltzmann weight factor for the energy
          
          if(bweight < R)
            {for(int i=0;i<N;i++)
                {trial[i]=lattice[i];}  // copying the old configuration into new trial configuration
                 continue;}
        
          else
            {for(int i=0;i<N;i++)
              {lattice[i]=trial[i];} } 
              
          if(t>500)  
            {
            minimumenergy+=E;
            flag++;
            }  
              }
 
      for(int i=0;i<N;i++)
        {trial[i]=lattice[i]; }   // copying the new configuration into new trial configuration
      

         
        
      data<<t<<" "<<Etrial<<"\n";    
  }  
      
  data.close();
  cout<<"\n"<<minimumenergy/flag;
   return 0; 
}