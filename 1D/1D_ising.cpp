#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <vector>
#include <cstdio>
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;

string folder_path = fs::current_path().string();

double energystrength=1;
int numberoftrials=10;
int N= 400 ;//Number of lattice sites 
int maxiter=10000;// number of runs
double Temperature=1; // Temperature
double B=1/Temperature ; // 1/KbT 
int flag=0;


vector<vector <double>> magnetization(numberoftrials),magnetization_time;
int seed =1278347;


random_device rd;
mt19937 gen(rd());

uniform_int_distribution choose_random_site(0,N-1);
uniform_real_distribution <double> boltzmannweight(0,1);





double energy( vector <int>& lattice)
{   double energy=0;      
    for ( int i=0 ; i<N ; i++ )
    {   if(i==N-1)
          {energy+= double(lattice[i]*lattice[0]);
          }
        else
          {energy+= double(lattice[i]*lattice[i+1]);
          }  
       
        if(i==0)
        {energy+= double(lattice[i]*lattice[N-1]);
          }
        else
          {energy+= double(lattice[i]*lattice[i+1]);
          }  
          
    } 
    
    energy*=-energystrength;  
    return energy;
}
void initialize_hot(vector <int>& lattice,vector <int>& trial_lattice)
{   srand(seed);
    for(int i=0;i<N;i++)
    {    
        if(rand()%10<=5)
          {lattice[i]=1;
           trial_lattice[i]=1;}
        else
          {lattice[i]=-1;
           trial_lattice[i]=-1;}
    }
    seed+=10;
}

void initialize_cold(vector <int>& lattice,vector <int>& trial_lattice)
{  
    for(int i=0;i<N;i++)
          {lattice[i]=-1;
           trial_lattice[i]=-1;}

    
}

double calculate_magnetization(vector <int>& lattice  )
{
  // calculating maghnetizzation vs time graph
    double magnet=0.0;
    for(int i=0;i<N;i++)  
        magnet+=lattice[i];
    magnet/=N;
  return magnet;
  
}

void savemagnetization(vector <double> magnetization,vector<int> time)
{

}



time_t trial_time,start_time=time(NULL) , finish_time;

int main(){
  for(int trial=0;trial<numberoftrials;trial++){
    trial_time=time(NULL);
    vector <int> lattice(N,0),trial_lattice(N,0);
    initialize_hot(lattice,trial_lattice);
    ofstream data(folder_path+"/data/systemimage_"+to_string(trial)+"_.txt");
    for (int t=0;t<maxiter;t++){  
      double E = 0.0;
      double Etrial=0.0;       
      
      E = energy(lattice);                       
      int choose=choose_random_site(gen);  // Choosing a random particle
      trial_lattice[choose]*=-1;                   // flipping the spin of the random particle
      Etrial=energy(trial_lattice);                // Calculating new configuration energy

      if(Etrial < E) {
        lattice=trial_lattice;       // checking if the trial energy is less than the original energy and if true accepting the new configuration
        for(int i=0;i<N;i++){data<<trial_lattice[i]<<" ";} data<<"\n"; 
        calculate_magnetization(lattice);         // calculating magetization
      }     
            
      else  // If the trial energy is greater then 
      { double R=boltzmannweight(gen);  // choosing a number between (0,1) 
        double bweight = exp(-(Etrial-E)*B); // calculating the boltzmann weight factor for the energy
        
        if(bweight >= R){ 
          lattice=trial_lattice;            // accepting the new configuration and storing the data 
          for(int i=0;i<N;i++){data<<trial_lattice[i]<<" ";} data<<"\n";    
          calculate_magnetization(lattice);         // calculating magetization
        }
      }
      trial_lattice=lattice;   // copying the new configuration into new trial configuration
      

    }  // end of time loop
    data.close();
    cout<<"\n"<<"Time to calculate trial = "<<trial_time-time(NULL)<<" seconds";
    
  }
      
}