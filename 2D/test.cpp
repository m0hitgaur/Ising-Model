#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <vector>
#include <cstdio>
#include <filesystem>
using namespace std;
namespace fs = filesystem;

string folder_path = fs::current_path().string();

double J=1;
int numberoftrials=10;
int Lx= 32,Ly= 32 ;//Number of lattice sites 
int maxiter=50000;// number of runs
double Tc=2.27;
double Temperature=0.5*Tc; // Temperature




int seed =1278347;


random_device rd;
mt19937 gen(rd());
uniform_real_distribution <double> boltzmannweight(0,1);
time_t trial_time,start_time=time(NULL) , finish_time;


double energy( vector <vector <int>>& lattice,int i,int j)
{   double energy=0;      
    energy+= lattice[i][(Ly+j+1)%Ly] + lattice[i][(Ly+j-1)%Ly] + lattice[(Lx+i+1)%Lx][j] + lattice[(Lx+i-1)%Lx][j];
    energy*=2*J*lattice[i][j];  
    return energy;
}

void initialize_hot(vector <vector <int>>& lattice)
{   srand(seed);
    for(int i=0;i<Lx;i++)
    { for(int j=0;j<Ly;j++){
        if(rand()%10<5)lattice[i][j]=1;
        else lattice[i][j]=-1;
      }
    }
    seed+=10;
}


void initialize_cold(vector <vector <int>>& lattice){  
  for(int i=0;i<Lx;i++)for(int j=0;j<Ly;j++)lattice[i][j]=-1;      
}


void magnetization_func(vector <vector <int>>& lattice,int t,int t_ ,int trial ,vector<vector<double>>& magnetization,vector<double>& magnetization_time ){
  // calculating maghnetizzation vs time graph  
    double magnet=0.0;
    for(const auto row:lattice)for(int element:row )  magnet+=element;
    magnet/=Lx*Ly;
    magnetization[trial][t_]=magnet;
    if(trial==0)magnetization_time[t_]=t;  
  }


double calculate_mag(vector <vector <int>>& lattice){
  double magnet=0.0;
  for(const auto row:lattice)for(int element:row )  magnet+=element;
  magnet/=Lx*Ly;
  return magnet;
}

void calculate_moments(vector <vector <int>>& lattice,int t_,int trial,vector <vector <double>>&secondmoment,vector <vector <double>>& fourthmoment)
{
  // calculating maghnetizzation vs time graph
  double phi=calculate_mag(lattice);
  secondmoment[trial][t_]=phi*phi;  
  fourthmoment[trial][t_]=phi*phi*phi*phi; 

}


void visualize(vector <vector <int>>& lattice,int t,int trial){
 fs::path data_dir_path = fs::path(folder_path) / "data";

// 2. Check if the directory exists
if (!fs::exists(data_dir_path)) {
    try {
        fs::create_directories(data_dir_path);
    } catch (const exception& e) {
        cerr << "Error creating directory " << data_dir_path << ": " << e.what() << endl;
    }
}
  ofstream file(folder_path+"/data/lattice_"+to_string(trial)+"_"+to_string(t)+"_.txt");
  for(auto row:lattice){for(int element:row) file<<element<<" "; file<<"\n";}
  file.close();
}

void calculate_U4(vector<double>& magnetization_time,vector <vector <double>>&secondmoment,vector <vector <double>>& fourthmoment){
  ofstream f(folder_path+"/U4_vs_time.txt");
  for(int t=0;t<size(magnetization_time);t++){
    double secondmom=0.0,fourthmom=0.0;
    for(int trial=0;trial<numberoftrials;trial++){
      secondmom+=secondmoment[trial][t];
      fourthmom+=fourthmoment[trial][t];
    }
    secondmom/=numberoftrials;
    fourthmom/=numberoftrials;
    double U4=1-(fourthmom/(3*secondmom*secondmom));
    f<<magnetization_time[t]<<" "<<U4<<"\n";
  }
  
  f.close();
}

void savemagnetization(vector<vector<double>>& magnetization,vector<double>& magnetization_time ){
    ofstream fil(folder_path+"/mag_vs_time.txt");
    vector<double> magnetization_avg(size(magnetization[0]),0);
    for(int i=0;i<size(magnetization[0]);i++){
        for(int tr=0;tr<numberoftrials;tr++) magnetization_avg[i]+=magnetization[tr][i];
        magnetization_avg[i]/=numberoftrials;
        fil<<magnetization_time[i]<<" "<<magnetization_avg[i]<<"\n";
    }
    fil.close();
}

int main(){
    int numberoftimesteps=0;
    for (int t=0;t<maxiter;t++){ 
        int tf;
        if (t<10) tf=1;
        if(t>10)  tf=10;
        if(t>100)  tf=50;
        if (t>1000)  tf=100;
        if (t>10000)  tf=500;
        if(t%tf==0){numberoftimesteps+=1;}
    }
    
    vector<vector<double>> magnetization(numberoftrials,vector<double>(numberoftimesteps,0));
    vector<double> magnetization_time(numberoftimesteps);
    vector<vector<double>> secondmoment(numberoftrials,vector<double>(numberoftimesteps,0));
    vector<vector<double>> fourthmoment(numberoftrials,vector<double>(numberoftimesteps,0));
    
  for(int trial=0;trial<numberoftrials;trial++){
    trial_time=time(NULL);
    vector<vector <int>> lattice(Lx,vector<int>(Ly,0));
    initialize_hot(lattice);
    int t_=0;
    for (int t=0;t<maxiter;t++){ 
      int tf;
      if (t<10) tf=1;
      if(t>10)  tf=10;
      if(t>100)  tf=50;
      if (t>1000)  tf=100;
      if (t>10000)  tf=500;
      if(t%tf==0){magnetization_func(lattice,t,t_,trial,magnetization,magnetization_time);calculate_moments(lattice,t_,trial,secondmoment,fourthmoment);t_++;} 
      if(t%10000==0)cout<<endl<<"<<"<<t;
      for(int n=0;n<Lx*Ly;n++){                                
        int choose_x=rand()%Lx,choose_y=rand()%Ly;  // Choosing a random site
        double dE=energy(lattice,choose_x,choose_y);                // Calculating new configuration energy

        // checking if the test energy is less than the original energy and if true accepting the new configuration and  calculating magetization
        if(dE<0){lattice[choose_x][choose_y]*=-1;} 
        // If the test energy is greater then then choose a random number between (0,1) 
        //and compare it with the boltzmann weight factor and accept the configuration if the weight factor is greater than the random number 
        else if(exp(-(dE)/Temperature)>=double(rand()%100)/100.0  )lattice[choose_x][choose_y]*=-1;
        
      }     

        
    }  
    cout<<"\n"<<"Time to calculate trial = "<<time(NULL)-trial_time<<" seconds";
    
  }

  savemagnetization(magnetization,magnetization_time);
  
  calculate_U4(magnetization_time,secondmoment,fourthmoment);
}
