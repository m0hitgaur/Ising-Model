#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include<filesystem>
#include <numeric> // Required for accumulate
namespace fs = std::filesystem;
using namespace std;


class IsingModel1D {
private:
    int N;                          // Number of spins
    vector <int> spins;         // Spin configuration (-1 or +1)
    double J;                       // Coupling ant (>0 for ferromagnetic)
    double T;                       // Temperature
    double beta;                    // 1/kT (k_B = 1)
    mt19937 gen;               // Random number generator
    uniform_real_distribution<> uniform_dist;
    uniform_int_distribution<> site_dist;
    
public:
    IsingModel1D(int size, double coupling, double temperature, unsigned int seed) 
        : N(size), J(coupling), T(temperature), beta(1.0/temperature),
          gen(seed), uniform_dist(0.0, 1.0), site_dist(0, size-1) {
        
        spins.resize(N);
        // Initialize with random spins
        for (int i = 0; i < N; i++) spins[i] = (uniform_dist(gen) < 0.5) ? -1 : 1;
    }
    
    // Initialize all spins up
    void initializeAllUp() {
        for (int i = 0; i < N; i++) {
            spins[i] = 1;
        }
    }
    
    // Initialize all spins down
    void initializeAllDown() {
        for (int i = 0; i < N; i++) {
            spins[i] = -1;
        }
    }
    
    // Calculate energy of the system
    double calculateEnergy()  {
        double energy = 0.0;
        for (int i = 0; i < N; i++) {
            // Periodic boundary conditions
            int right = (i + 1) % N;
            energy -= J * spins[i] * spins[right];
        }
        return energy;
    }
    
    // Calculate magnetization of the system
    double calculateMagnetization()  {
        double mag = 0.0;
        for (int i = 0; i < N; i++) {
            mag += spins[i];
        }
        return mag / N;  // Normalized magnetization
    }
    
    // Calculate local energy for spin at site i
    double localEnergy(int i)  {
        int left = (i - 1 + N) % N;
        int right = (i + 1) % N;
        return -J * spins[i] * (spins[left] + spins[right]);
    }
    
    // Metropolis algorithm for single spin flip
    void metropolisStep() {
        // Choose random site
        int i = site_dist(gen);
        
        // Calculate energy change if we flip this spin
        double E_old = localEnergy(i);
        double dE = -2 *E_old;
        
        // Metropolis acceptance criterion
        if (dE <= 0 || uniform_dist(gen) < exp(-beta * dE)) {
            spins[i] *= -1;  // Flip the spin
        }
    }
    
    // Perform one Monte Carlo sweep (N attempted flips)
    void monteCarloSweep() {
        for (int i = 0; i < N; i++) {
            metropolisStep();
        }
    }
    
    // Get spin configuration
    vector<int> getSpins()  {
        return spins;
    }


    // Run simulation for given number of sweeps
    void simulate( int totalsweeps,double& avg_energy, double& avg_mag, double& avg_mag_abs,double& heat_capacity, double& susceptibility, int& run) {
        
        double sum_E = 0.0, sum_E2 = 0.0;
        double sum_M = 0.0, sum_M2 = 0.0;
        vector <double> mag_vs_time;
        vector <int> time_mag;
        int sweep_flag=0;
        for (int sweep = 0; sweep < totalsweeps; sweep++) {
            ofstream file("/data/run_"+to_string(run)+"/lattice_"+to_string(sweep)+".dat");
            monteCarloSweep();
            double E = calculateEnergy();
            double M = calculateMagnetization() * N; 
            sum_E += E;
            sum_E2 += E * E;
            sum_M += M;
            sum_M2 += M * M;
            
            if(sweep<10)sweep_flag=1;
            else if(sweep>10)sweep_flag=10;
            else if(sweep>100)sweep_flag=100;
            else if(sweep>1000)sweep_flag=500;
            
            if(sweep%sweep_flag==0){
                mag_vs_time.push_back(M);
                time_mag.push_back(sweep);
                vector<int> lattice=getSpins();
                for(auto a:lattice)file<<a<<" ";
            }
            file.close();
        }
        
        // Calculate averages
        double E_avg = sum_E / totalsweeps;
        double M_avg = sum_M / totalsweeps;

        avg_energy = E_avg / N;
        avg_mag = M_avg / N;
                
        double E2_avg = sum_E2 / totalsweeps;
        double M2_avg = sum_M2 / totalsweeps;
        
        
        heat_capacity = (E2_avg - E_avg * E_avg) / (T * T * N);
        susceptibility = (M2_avg - M_avg * M_avg) / (T * N);
    }
    
    // Set temperature
    void setTemperature(double temp) {
        T = temp;
        beta = 1.0 / temp;
    }
    
    // Print current configuration
    void printConfiguration()  {
        for (int i = 0; i < N; i++) {
            cout << (spins[i] > 0 ? "+" : "-") << " ";
        }
        cout << endl;
    }
    

};


double calculate_mean( vector<double>& v) {
    if (v.empty()) return 0.0;
    double sum = accumulate(v.begin(), v.end(), 0.0);
    return sum / v.size();
}

double calculate_std_dev( vector<double>& v) {
    if (v.size() < 2) return 0.0;
    double mean = calculate_mean(v);
    double sq_sum = 0.0;
    for ( auto& x : v) {
        sq_sum += (x - mean) * (x - mean);
    }
    return sqrt(sq_sum / (v.size()));
}


void make_dir(int& num_runs){
    fs::path folder_path = fs::current_path();
    cout<<"CREATING OUTPUT DIRECTORIES";
    vector<string>directories;
    directories.push_back(folder_path.string()+"/data");
    for(int run=0;run<num_runs;run++)directories.push_back(folder_path.string()+"/data/run_"+to_string(run));

    for(auto dir_path:directories){
        cout<<dir_path<<"\n";
        if (fs::exists(dir_path)) {
            // The path exists, now check if it's a directory
            if (fs::is_directory(dir_path)) {
                //std::cout << "Success! The directory already exists." << std::endl;
            } else {
                // This handles the case where a file with the same name exists
                //std::cout << "Error: A file with the same name already exists at this path." << std::endl;
            }
        } else {
            // The directory does not exist, so we create it
            //std::cout << "Directory does not exist. Creating it..." << std::endl;
            try {
                // Attempt to create the directory
                if (fs::create_directory(dir_path)) {
                    //std::cout << "Directory successfully created!" << std::endl;
                } else {
                    cout << "Failed to create the directory." << std::endl;
                }
            } catch (const fs::filesystem_error& e) {
                // Catch any errors during directory creation (e.g., permissions issues)
                cerr << "Filesystem error: " << e.what() << std::endl;
            }
        }
    }
    cout<<"OUTPUT DIRECTORIES CREATED SUCCESSFULLY \n";

}

int main() {
    // System parameters
     int N = 100;                // Number of spins
     double J = 1.0;             // Coupling const
     int numberofsweeps = 5000;   // total sweeps 
     int num_runs = 10;          // Number of independent runs to average over

     make_dir(num_runs);
    // Use a random device to seed the random number generator for each run
    random_device rd;

    // --- Main Simulation ---
    cout << "1D Ising Model Simulation\n";
    cout << "System size: " << N << " spins\n";
    cout << "Coupling: J = " << J << "\n";
    cout << "Averaging over " << num_runs << " independent runs per temperature.\n\n";
    
    cout << fixed << setprecision(5);
    cout << setw(8) << "T"
         << setw(14) << "Energy/spin"
         << setw(14) << "|M|"
         << setw(14) << "Heat Cap"
         << setw(14) << "Suscept" << endl ;   
    cout << string(64, '-') << endl;
    
    // Output file for plotting
    ofstream outfile("ising_1d_averaged_results.txt");
    outfile << "# T Energy E_err Mag_abs M_err HeatCap C_err Suscept Chi_err\n";
    
    // Temperature range
    vector<double> temperatures;
    for (double T = 0.1; T <= 1; T += 0.1) {
        temperatures.push_back(T);
    }
    
    for (double T : temperatures) {
        // Store results from each run for the current temperature
        vector<double> run_energies;
        vector<double> run_magnitudes;
        vector<double> run_heat_caps;
        vector<double> run_suscepts;

        // Perform multiple independent runs
        for (int run = 0; run < num_runs; run++) {
            // Create a new model with a unique seed for each run
            IsingModel1D model(N, J, T, rd());
            
            double avg_energy, avg_mag, avg_mag_abs, heat_cap, suscept;
            model.simulate(numberofsweeps ,avg_energy, avg_mag, avg_mag_abs, heat_cap, suscept,run);
            // Store the results of this run
            run_energies.push_back(avg_energy);
            run_magnitudes.push_back(avg_mag_abs);
            run_heat_caps.push_back(heat_cap);
            run_suscepts.push_back(suscept);
        }

        // Calculate final averages and standard errors over the runs
        double final_avg_energy = calculate_mean(run_energies);
        double stderr_energy = calculate_std_dev(run_energies) / sqrt(num_runs);

        double final_avg_mag = calculate_mean(run_magnitudes);
        double stderr_mag = calculate_std_dev(run_magnitudes) / sqrt(num_runs);

        double final_avg_heat_cap = calculate_mean(run_heat_caps);
        double stderr_heat_cap = calculate_std_dev(run_heat_caps) / sqrt(num_runs);

        double final_avg_suscept = calculate_mean(run_suscepts);
        double stderr_suscept = calculate_std_dev(run_suscepts) / sqrt(num_runs);
        
        cout << setw(8) << T
                  << setw(14) << final_avg_energy
                  << setw(14) << final_avg_mag
                  << setw(14) << final_avg_heat_cap
                  << setw(14) << final_avg_suscept << endl;
        
        outfile << T << " " << final_avg_energy << " " << stderr_energy
                << " " << final_avg_mag << " " << stderr_mag
                << " " << final_avg_heat_cap << " " << stderr_heat_cap
                << " " << final_avg_suscept << " " << stderr_suscept << "\n";
    }
    
    outfile.close();
    
    cout << "\nResults saved to 'ising_1d_averaged_results.txt'\n"; 

    
    return 0;
}