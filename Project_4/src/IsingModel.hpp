#ifndef ISINGMODEL
#define ISINGMODEL

#include <armadillo>

class IsingModel{
    private: 
        arma::mat m;
        int L;
        arma::vec delta_Energy_arr;
        int numb_accepted_states;
        int numb_spins;
        double temperature;
        double energy;
        double magnetism; 
        
        //Initialize them here so we can use them in the whole class 
        double energy_mean;
        double energy2_mean;
        double magnetism_mean;
        double magnetism2_mean;
        double magnetism_abs_mean;

        double specific_heat;
        double susceptibility;


    public:
        
        
        arma::vec epsilon_arr;
        
        IsingModel(int latice_dimentions);
        void initialise(double temperature_in);
        void initialise(double temperature_in, arma::mat m_in);
        double totalEnergy();
        double totalMagnetization();
        int idx(int i);
        double deltaEnergy(int i, int j);

        //Methods for markov chain approach
        void metropolis();
        void metropolisOMP(arma::mat &matrix, double &E, double &M, arma::vec dE_arr, int &numb_accepted);
        void monteCarlo(int cycles);
        void monteCarloBurn(int cycles, arma::vec &num_cycles, arma::vec &energy_vec, arma::vec &m_vec);
        void monteCarloOMP(int cycles);



        double deltaEnergy(int i, int j, arma::mat matrix);


        void burnInTime();

        //Outprint
        void write_results();
        double get_energy_mean();
        double get_energy2_mean();
        double get_m_abs_mean();
        double get_m2_abs_mean();
        double get_numb_accepted_states();
        double get_heat(); 
        double get_sus(); 
        void generate_error(double &E_mean, double &E2_mean, double &m_abs_mean, double &m2_mean, double &spec_heat, double &susc);
        
};
#endif