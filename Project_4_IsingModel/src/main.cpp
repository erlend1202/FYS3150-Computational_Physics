#include <iostream>
#include <armadillo>
#include "IsingModel.hpp"
#include <omp.h>

#include <chrono>

#include "utilities.hpp"
#include "omp_rng.hpp"

void totalEnergyTest1();
void test_montecarlo();
void test_montecarlo_OMP();
void approx_prob();
void compare_burn_in_time();
void time_test_OMP();
void laticeSizeComp();
void error_analysis();
void generate_energy_and_magnetism();
void test_new_rand();

int main(){
    int seed = 123123;
    omp_rng::initialize_omp_rng_container(seed);  // Initialization of random number generator
    compare_burn_in_time();
    return 0;
}


/**
 * Method for testing speedup of the paralellization 
 * of the Ising Model
 * Testing both with Latice size 20 and 500000 cycles of montecarlo
**/
void time_test_OMP(){
    int L = 20;
   
    int cycles = 500000;

    auto t1 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < 8; i ++){
        std::cout << i << std::endl;
        IsingModel test = IsingModel(L);
        arma::mat m_in = makeMatrixRandom(L);
        test.initialise(1.0, m_in);
        test.monteCarlo(cycles); 
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    double duration_seconds_normal = std::chrono::duration<double>(t2 - t1).count()/10;


    t1 =  std::chrono::high_resolution_clock::now();
    #pragma omp parallel for 
    for (size_t i = 0; i < 8; i ++){
        std::cout << i << std::endl;
        IsingModel test = IsingModel(L);
        arma::mat m_in = makeMatrixRandom(L);
        test.initialise(1.0, m_in);
        test.monteCarlo(cycles); 
    }
    t2 = std::chrono::high_resolution_clock::now();
    double duration_seconds_OMP = std::chrono::duration<double>(t2 - t1).count()/10;

    std::cout << "time normal: " << duration_seconds_normal << "  time OMP: " 
    << duration_seconds_OMP << std::endl;
}

/**
 * Function that generates the expected energy and magnetization per spin for different numbers of 
 * cycles with the monte carlo simulation and different temperatures T = 1 and T = 2.4
**/
void compare_burn_in_time(){
    int L = 20;
    double T1 = 1;
    double T2 = 2.4;

    IsingModel rand_model1 = IsingModel(L);
    IsingModel rand_model2 = IsingModel(L);

    IsingModel ord_model1 = IsingModel(L);
    IsingModel ord_model2 = IsingModel(L);

    arma::mat rand_spins = makeMatrixRandom(L);
    arma::mat ord_spins = makeMatrixOrdered(L);

    rand_model1.initialise(T1, rand_spins);
    rand_model2.initialise(T2, rand_spins);

    ord_model1.initialise(T1, ord_spins);
    ord_model2.initialise(T2, ord_spins);

    int n = 1000000;
    int start = 1; 

    arma::vec rand_energy1 = arma::vec(n);
    arma::vec rand_energy2 = arma::vec(n);

    arma::vec ord_energy1 = arma::vec(n);
    arma::vec ord_energy2 = arma::vec(n);

    arma::vec rand_abs_m1 = arma::vec(n);
    arma::vec rand_abs_m2 = arma::vec(n);

    arma::vec ord_abs_m1 = arma::vec(n);
    arma::vec ord_abs_m2 = arma::vec(n);

    arma::vec cycles = arma::vec(n);

    rand_model1.monteCarloBurn(n, cycles, rand_energy1, rand_abs_m1);

    rand_model2.monteCarloBurn(n, cycles, rand_energy2, rand_abs_m2);

    ord_model1.monteCarloBurn(n, cycles, ord_energy1, ord_abs_m1);

    ord_model2.monteCarloBurn(n, cycles, ord_energy2, ord_abs_m2);

    
    print_3_vecs(rand_energy1, rand_abs_m1, cycles, n, "T1_mean_energy_and_abs_m_random.txt", 15);
    print_3_vecs(ord_energy1, ord_abs_m1, cycles, n, "T1_mean_energy_and_abs_m_ordered.txt", 15);

    print_3_vecs(rand_energy2, rand_abs_m2, cycles, n, "T2_mean_energy_and_abs_m_random.txt", 15);
    print_3_vecs(ord_energy2, ord_abs_m2, cycles, n, "T2_mean_energy_and_abs_m_ordered.txt", 15);

}

/**
 * Writes the energy and magnetic information from 1000000
 * different cycles for monteCarloBrun to the file error_analysis.txt
 */
void error_analysis(){
    int L = 2;
    double T = 1;

    IsingModel model = IsingModel(L);
    arma::mat m_in = makeMatrixRandom(L);
    model.initialise(T, m_in);
    std::cout << 1 << std::endl;

    arma::mat spins = makeMatrixRandom(L);
    std::cout << 2 << std::endl;

    int n = 1000000;

    arma::vec energy_vec = arma::vec(n);
    std::cout << 3 << std::endl;

    arma::vec abs_m_vec = arma::vec(n);
    std::cout << 4 << std::endl;

    arma::vec cycles_vec = arma::vec(n);
    std::cout << 5 << std::endl;

    model.monteCarloBurn(n, cycles_vec, energy_vec, abs_m_vec);
    std::cout << 6 << std::endl;

    print_3_vecs(energy_vec, abs_m_vec, cycles_vec, n, "error_analysis.txt", 15);
    std::cout << 7 << std::endl;
}

/**
 * Writes energy, magnetism, susceptibility and specific heat data to files
 * for 10000 cycles of montecarlo.
 */
void generate_energy_and_magnetism(){
    int L = 2;
    int cycles = 100000;
    int N = 40;
    
    arma::vec E_mean_vec = arma::vec(N);
    arma::vec E2_mean_vec = arma::vec(N);
    arma::vec m_abs_mean_vec = arma::vec(N);
    arma::vec m2_mean_vec = arma::vec(N);
    arma::vec spec_heat_vec = arma::vec(N);
    arma::vec susc_vec = arma::vec(N);

    arma::vec temp = arma::vec(N);

    double E_mean, E2_mean, m_abs_mean, m2_mean, spec_heat, susc;

    double stepsize = (double)3/N; 
    double tmp = 1; 
    for (size_t i = 0; i < N; i ++){
        IsingModel model = IsingModel(L);
        arma::mat m_in = makeMatrixRandom(L);
        model.initialise(tmp, m_in);
        model.monteCarlo(cycles);
        model.generate_error(E_mean, E2_mean, m_abs_mean, m2_mean, spec_heat, susc);

        E_mean_vec(i) = E_mean;
        E2_mean_vec(i) = E2_mean;

        m_abs_mean_vec(i) = m_abs_mean;
        m2_mean_vec(i) = m2_mean;

        spec_heat_vec(i) = spec_heat;
        susc_vec(i) = susc;

        temp(i) = tmp; 
        tmp += stepsize; 
    }

    print_3_vecs(E_mean_vec, E2_mean_vec, temp, N, "numerical_energy_T.txt", 15);
    print_3_vecs(m_abs_mean_vec, m2_mean_vec, temp, N, "numerical_magnetism_T.txt", 15);
    print_3_vecs(spec_heat_vec, susc_vec, temp, N, "numerical_spec_susc_T.txt", 15);

    std::cout << temp << std::endl; 
}
/**
 * Function for testing our monte carlo implementation against
 * the analytical solutions we found.
**/
void test_montecarlo(){
    int L = 2;
    IsingModel test = IsingModel(L);
    arma::mat m_in = makeMatrixRandom(L);

    test.initialise(1.0, m_in);
    for (int i = 10; i <= 1000000; i*=10){
        std::cout << "cycles = " << i << std::endl;
        test.monteCarlo(i); 

        //function to write the results from our monte carlo problem
        test.write_results();
    }
}


//for problem 6

/**
 * Generating values for two different temeperatures 
 * and saves the results to the file problem6_test.txt
 */
void approx_prob(){
    int L = 20;
    arma::vec T_1;
    arma::vec T_24;
    IsingModel test = IsingModel(L);

    int cycles = 1000;
    //simulating for T=1.0
    test.initialise(1.0);
    test.burnInTime();
    test.monteCarlo(cycles);
    T_1 = test.epsilon_arr;

    std::cout << "halfways" << std::endl;

    //simulating for T=2.4
    test.initialise(2.4);
    test.burnInTime();
    test.monteCarlo(cycles);
    T_24 = test.epsilon_arr;

    std::cout << "done" << std::endl;

    print_2_vecs(T_1, T_24, cycles, "problem6_test.txt", 15);
}

/**
 * Function for testing energy, magnetization and delta energy calculations
 * made by the Ising model class
 */
void totalEnergyTest1(){
    int L = 5;
    arma::mat b(L,L);
    b.fill(1);
    b(1,1) = -1;
    b(1,0) = -1;
    b(1,2) = -1;

    std::cout << b << std::endl;

    IsingModel test = IsingModel(L);
    test.initialise(1.0, b);
    std::cout << "total energi: " << test.totalEnergy() << std::endl;
    std::cout << "total magnetization: " << test.totalMagnetization() << std::endl;
    std::cout << "deltaE: " << test.deltaEnergy(1,1) << std::endl;
}

/**
 * Function calculating and extracting data from Ising models
 * of different sizes
 * Cycles = 2000000
 * latsizes = 40, 60, 80, 100
 */
void laticeSizeComp(){
    arma::vec latSize = arma::vec(4);
    latSize(0) = 100;
    latSize(1) = 80;
    latSize(2) = 60;
    latSize(3) = 40;



    int num_cycles = 2000000;

    arma::vec temp1 = deltaSpace(2.1,2.25,0.01);
    arma::vec temp2 = deltaSpace(2.255,2.38,0.005);
    arma::vec temp3 = deltaSpace(2.39,2.401,0.01);

    arma::vec temp = arma::join_cols(temp1, temp2, temp3);
    std::cout << temp << std::endl;
    auto t1 = std::chrono::high_resolution_clock::now();

    
    for (int i = 0; i < latSize.size(); i++){
        std::cout << "latsize" << latSize(i) << std::endl;
        std::vector<arma::vec> data;


        arma::vec eps = arma::vec(temp.size());
        arma::vec eps_2 = arma::vec(temp.size());
        arma::vec mag = arma::vec(temp.size());
        arma::vec mag_2 = arma::vec(temp.size());
        arma::vec heat = arma::vec(temp.size());
        arma::vec susc = arma::vec(temp.size());
        arma::vec temperature = arma::vec(temp.size());
        arma::vec numAccStates = arma::vec(temp.size());
        
        #pragma omp parallel
        {   
            #pragma omp for
            for (int j = 0; j < temp.size(); j++){
                std::cout << j << "/" << temp.size() << std::endl;
                IsingModel model(latSize[i]);
                model.initialise(temp[j]);
                model.burnInTime();
                model.monteCarlo(num_cycles);                
                
                eps[j] = model.get_energy_mean()/(latSize(i)*latSize(i));
                eps_2[j] = model.get_energy2_mean();

                mag[j] = model.get_m_abs_mean();
                mag_2[j] = model.get_m2_abs_mean();

                heat[j] = model.get_heat();
                susc[j] = model.get_sus();

                temperature[j] = temp[j];
                numAccStates[j] = model.get_numb_accepted_states();
            }
            #pragma omp barrier
        }

        data.push_back(eps);
        data.push_back(mag);
        data.push_back(heat);
        data.push_back(susc);
        data.push_back(temperature);
        data.push_back(eps_2);
        data.push_back(mag_2);
        data.push_back(numAccStates);

        std::string fileName = "Mat" + std::to_string((int)latSize(i)) + "info.txt";
        print_n_vecs(data, temp.size(), fileName.c_str(), 15);
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    double duration_seconds = std::chrono::duration<double>(t2 - t1).count();
    std::cout << "Time = " << duration_seconds << std::endl;
}

void test_new_rand(){
    std::cout << omp_rng::get_random_int_0_n(300) << std::endl;
    std::cout << omp_rng::get_random_int_0_n(300) << std::endl;
    std::cout << omp_rng::get_random_int_0_n(10) << std::endl;
    omp_rng::initialize_omp_rng_container(1);
    std::cout << omp_rng::get_random_int_0_n(100) << std::endl;
    omp_rng::initialize_omp_rng_container(1);
    std::cout << omp_rng::get_random_int_0_n(1000) << std::endl;

}