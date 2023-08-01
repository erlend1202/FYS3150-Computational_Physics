#include <iostream>
#include <armadillo>
#include "IsingModel.hpp"
#include <cmath>
#include <omp.h>
#include "utilities.hpp"
#include "omp_rng.hpp"


#include <omp.h>


/**
 * Initializer for the class IsingModel.
 * ----------
 * Parameters
 * ----------
 * int latice_dimensions: size of our spin-matrix dimensions
**/
IsingModel::IsingModel(int latice_dimensions){
    L = latice_dimensions;
    numb_spins = L*L;
    energy = 0;
    magnetism = 0;
}

/**
 * Method to start our system and set a temperature T.
 * ----------
 * Parameters
 * ----------
 * double temperature_in: temperature we want to set for our system
**/
void IsingModel::initialise(double temperature_in){
    temperature = temperature_in;
    numb_accepted_states = 0;  

    delta_Energy_arr.set_size(17);
    delta_Energy_arr.fill(0);

    for (int i = -8; i <= 8; i += 4){
        delta_Energy_arr[i + 8] = exp(-i/temperature);
    }

    srand(time(NULL));
    m = makeMatrixRandom(L);

    energy = totalEnergy();
    magnetism = totalMagnetization();

}

//overloaded to take in a matrix of choice
void IsingModel::initialise(double temperature_in, arma::mat m_in){
    temperature = temperature_in;
    numb_accepted_states = 0;  

    delta_Energy_arr.set_size(17);
    delta_Energy_arr.fill(0);

    for (int i = -8; i <= 8; i += 4){
        delta_Energy_arr[i + 8] = exp(-i/temperature);
    }

    srand(time(NULL));
    m = m_in;

    energy = totalEnergy();
    magnetism = totalMagnetization();
}

/**
 * Method for finding indexes using
 * periodic boundary conditions if index is out of range.
 * Returns the calculated index. 
 * ----------
 * Parameters
 * ----------
 * int i: value of our x or y index.
**/
int IsingModel::idx(int i){
    return (i + L) % L;
}


/**
 * Method for calculating the systems total energy.
 * Returns the total energy.
**/
double IsingModel::totalEnergy(){
    double sum = 0;
    //Center cases
    for (int i = 0; i < L-1; i++){
        for (int j = 0; j < L-1; j++){
            sum += m(i, j) * m(i+1, j);
            sum += m(i, j) * m(i, j+1);
        }
    }
    //Edge cases
    for (int k = 0; k < L-1; k++){
        //periodic boundary conditions
        sum += m(k,0)*m(k,L-1);
        sum += m(0,k)*m(L-1,k);

        //rightside pairs we missed
        sum += m(L-1,k)*m(L-1,k+1);

        //downward pairs we missed
        sum += m(k,L-1)*m(k+1,L-1);
    }

    //Periodic boundary condition on last element
    sum += m(L-1,L-1)*m(0,L-1);
    sum += m(L-1,L-1)*m(L-1,0);

    return -sum;
}



/**
 * Method for calculating the systems total magnetization.
 * Returns the total magnetization.
**/
double IsingModel::totalMagnetization(){
    double total = 0;
    for (int i = 0; i < m.n_cols; i++){
        for (int j = 0; j < m.n_rows; j++){
            total += m(i, j);
        }
    }
    return total;
}


/**
 * Calculates the change in energy using the index of the 
 * spin we configured.
 * This will calculate the change in energy that occurs IF we cinfigure
 * a spin, meaning you have to calculate this before flipping a spin,
 * otherwise you can use the negative of this to calculate the difference
 * in energy.
 * 
 * ----------
 * Parameters
 * ----------
 * int i: i-th index
 * int j: j-th index
 * 
 * Returns change in energy.
**/
double IsingModel::deltaEnergy(int i, int j){
    int up = m(idx(i+1), idx(j));
    int down = m(idx(i-1), idx(j));
    int right = m(idx(i), idx(j+1));
    int left = m(idx(i), idx(j-1));

    return 2*(up + down + right + left)*m(idx(i), idx(j));
}



//Markov Chain Monte Carlo approach methods.
/**
 * Method implementing the Markov Chain method using the
 * monte carlo approach.
 * ----------
 * Parameters
 * ----------
 * int cycles: amount of times we want to repeat the monte carlo simulation.
**/
void IsingModel::monteCarlo(int cycles){
    energy_mean = 0;
    energy2_mean = 0;
    magnetism_mean = 0;
    magnetism2_mean = 0;
    magnetism_abs_mean = 0;

    specific_heat = 0;
    susceptibility = 0;

    //for problem 6
    epsilon_arr.set_size(cycles);
    epsilon_arr.fill(0);


    for (int i = 0; i < cycles; i++){
        metropolis();

        energy_mean += energy;
        energy2_mean += energy*energy;
        magnetism_mean += magnetism;
        magnetism2_mean += magnetism*magnetism;
        magnetism_abs_mean += fabs(magnetism);

        // epsilon_arr(i) = energy/(double)numb_spins;
    }

    double normalize = double(1.0)/ ((double) cycles);
    energy_mean *= normalize;
    energy2_mean *= normalize;
    magnetism_mean *= normalize;
    magnetism2_mean *= normalize;
    magnetism_abs_mean *= normalize;
    

    double epsilon_mean = energy_mean/numb_spins;
    double epsilon2_mean = energy2_mean/numb_spins;
    double m_abs_mean = magnetism_abs_mean/numb_spins;
    double m_mean = magnetism_mean/numb_spins;
    double m2_mean = magnetism2_mean/numb_spins;


    specific_heat = (double)(1/(temperature*temperature)) * (energy2_mean - (energy_mean*energy_mean)) / (double)numb_spins;
    susceptibility = (double)(1/temperature) * (magnetism2_mean - (magnetism_abs_mean*magnetism_abs_mean)) / (double)numb_spins;
}

//Markov Chain Monte Carlo approach methods. Overloaded to save energy and magnetism for each cycle.
/**
 * Method implementing the Markov Chain method using the
 * monte carlo approach.
 * ----------
 * Parameters
 * ----------
 * int cycles: amount of times we want to repeat the monte carlo simulation.
**/
void IsingModel::monteCarloBurn(int cycles, arma::vec &num_cycles, arma::vec &energy_vec, arma::vec &m_vec){
    energy_mean = 0;
    energy2_mean = 0;
    magnetism_mean = 0;
    magnetism2_mean = 0;
    magnetism_abs_mean = 0;

    specific_heat = 0;
    susceptibility = 0;

    //for problem 6
    epsilon_arr.set_size(cycles);
    epsilon_arr.fill(0);

    for (int i = 0; i < cycles; i++){
        metropolis();

        energy_mean += energy;
        energy2_mean += energy*energy;
        magnetism_mean += magnetism;
        magnetism2_mean += magnetism*magnetism;
        magnetism_abs_mean += fabs(magnetism);

        energy_vec(i) = energy_mean*(double) 1.0 / ((double) i+1);
        m_vec(i) = magnetism_abs_mean*(double) 1.0 / ((double) i+1);
        num_cycles(i) = i; 

        epsilon_arr(i) = energy/(double)numb_spins;
    }

    double normalize = double(1.0)/ ((double) cycles);
    energy_mean *= normalize;
    energy2_mean *= normalize;
    magnetism_mean *= normalize;
    magnetism2_mean *= normalize;
    magnetism_abs_mean *= normalize;
    

    double epsilon_mean = energy_mean/numb_spins;
    double epsilon2_mean = energy2_mean/numb_spins;
    double m_abs_mean = magnetism_abs_mean/numb_spins;
    double m_mean = magnetism_mean/numb_spins;
    double m2_mean = magnetism2_mean/numb_spins;

    specific_heat = (double)(1/(temperature*temperature)) * (energy2_mean - (energy_mean*energy_mean)) / (double)numb_spins;
    susceptibility = (double)(1/temperature) * (magnetism2_mean - (magnetism_abs_mean*magnetism_abs_mean)) / (double)numb_spins;
}


/**
 * Method for returning generated values to be compared
 * to analytical answers.
**/
void IsingModel::generate_error(double &E_mean, double &E2_mean, double &m_abs_mean, double &m2_mean, double &spec_heat, double &susc){
    E_mean = energy_mean;
    E2_mean = energy2_mean; 
    m_abs_mean = magnetism_abs_mean;
    m2_mean = magnetism2_mean;
    spec_heat = specific_heat;
    susc = susceptibility;
}

/**
 * Method for writing out our results from the markov
 * chain monte carlo method.
**/
void IsingModel::write_results(){
    std::cout << "accepted states = " << numb_accepted_states << std::endl;
    std::cout << "expected energy: analytical = -7.936, numeric = " << energy_mean << std::endl;
    std::cout << "expected^2 energy: analytical = 63.87, numeric = "  << energy2_mean << std::endl;
    
    std::cout << "expected magnetism: analytical = 0, numeric = " << magnetism_mean << std::endl;
    std::cout << "expected abs magnetism: analytical = 3.99462, numeric = " << magnetism_abs_mean << std::endl;
    std::cout << "expected^2 magnetism: analytical = 15.973, numeric = " << magnetism2_mean << std::endl;

    std::cout << "specific_heat: analytical = 0.0321469, numeric = " << specific_heat << std::endl;
    std::cout << "susceptibility: analytical = 0.00401074, numeric = " << susceptibility << std::endl;

    std::cout << "" << std::endl;
}


/**
 * Method that implements and instance of the monte carlo approach
 * using the Metropolis-Hastings algorithm.
**/
void IsingModel::metropolis(){
    int x_idx;
    int y_idx;
    int dE;

    for (int i = 0; i < numb_spins; i++){

        x_idx = omp_rng::get_random_int_0_n(L-1);
        y_idx = omp_rng::get_random_int_0_n(L-1);

        dE = deltaEnergy(x_idx, y_idx);
        if (((double)omp_rng::get_random()) <= delta_Energy_arr[dE + 8]){

            m(x_idx, y_idx) *= -1;
            magnetism += (double) 2 * m(x_idx, y_idx);
            energy += (double) dE;
            numb_accepted_states++;
        }
    }
}

/**
 * Method for burning in our grid of spins.
 * Runs metropolis method 100000 times.
 */
void IsingModel::burnInTime(){
    for (int i = 0; i < 100000; i++){
        metropolis();
    }
}




double IsingModel::deltaEnergy(int i, int j, arma::mat matrix){
    int up = matrix(idx(i+1), idx(j));
    int down = matrix(idx(i-1), idx(j));
    int right = matrix(idx(i), idx(j+1));
    int left = matrix(idx(i), idx(j-1));

    return 2*(up + down + right + left)*matrix(idx(i), idx(j));
}



double IsingModel::get_energy_mean() {return energy_mean;}
double IsingModel::get_energy2_mean() {return energy2_mean;}
double IsingModel::get_m_abs_mean() {return magnetism_abs_mean / numb_spins;}
double IsingModel::get_m2_abs_mean() {return magnetism2_mean;}
double IsingModel::get_numb_accepted_states() {return numb_accepted_states;}
double IsingModel::get_heat() {return specific_heat;}
double IsingModel::get_sus() {return susceptibility;}

