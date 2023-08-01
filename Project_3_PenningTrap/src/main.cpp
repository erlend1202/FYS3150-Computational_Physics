#include <armadillo>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

#include "penningTrap.cpp"
#include "utilities.hpp"

void Particles_escaped();
void particleSimulationRandom(int seed, double dt, double totalTime, bool particleInteractions);
void particleSimulationSpecific(std::vector<Particle> particles, double dt, double totalTime, bool particleInteractions);

int main(){
    double q = 2;
    double m = 40.07;
    double d = pow(10, 4);

    //Partikkel 1
    /*arma::vec r1 = arma::vec(3).randn() * 0.1 * d;  // random initial position
    arma::vec v1 = arma::vec(3).randn() * 0.1 * d;  // random initial velocity


    //Partikkel 2
    arma::vec r2 = arma::vec(3).randn() * 0.1 * d;  // random initial position
    arma::vec v2 = arma::vec(3).randn() * 0.1 * d;  // random initial velocity

    Particle particle_1 = Particle(q, m, r1, v1);
    Particle particle_2 = Particle(q, m, r2, v2);

    std::vector<Particle> particles_in = {particle_1, particle_2};
    particleSimulationSpecific(particles_in, 0.01, 100, true);*/

    Particles_escaped();

    return 0;
}


/**
 * Initializes a penning trap with 100 particles for different amplitudes and frequencies.
 * Prints the amount of escaped particles over time to the files f_omega_coulomb_(amplitude).txt.
**/
void Particles_escaped(){

    //Problem 10
    double q = 2;
    double q1 = 2;
    double q2 = 2;

    double m = 40.07;

    double B_0 = 1;
    double V_0 = 0.0025;
    double d = 5*pow(10, 2);

    int num_particles = 100;


    // Set frequency
    double f_list[3] = {0.1, 0.4, 0.7};

    double omega_V;
    double f;

    arma::vec r;
    arma::vec v;

    // choose timestep and duration of simulation
    double delta = 0.1;
    int time_steps = 500/delta;
    int size_omega = 2.3/0.02;

    PenningTrap trapRK = PenningTrap(B_0, V_0, d);
    int index;

    arma::vec num_trapped = arma::vec(115);

    for (size_t f_it = 0; f_it < 3; f_it ++){
        f = f_list[f_it];
        index = 0;
        for (double omega_V = 0.2; omega_V <= 2.5; omega_V += 0.02){

            trapRK.particles.clear();

            for (size_t i = 0; i < num_particles; i++){
                r = arma::vec(3).randn() * 0.1 * d;  // random initial position
                v = arma::vec(3).randn() * 0.1 * d;  // random initial velocity

                Particle particle = Particle(q, m, r, v);
                trapRK.add_particle(particle);
            }
            //std::cout << trapRK.particles.size() << std::endl;

            double currentTime = 0;

            for (int i = 0; i < time_steps; i++){
                currentTime += delta;
                trapRK.evolve_RK4(delta, currentTime, f, omega_V);
            }

            num_trapped(index) = trapRK.num_escaped_particles();
            index += 1;
        }
        print_f_omega(f, omega_V, num_trapped, size_omega, time_steps, "f_omega_coulomb_" + std::to_string(f) + ".txt");
    }
}


/**
 * Simulates two randomly placed particles in the penning trap
 * @param seed The seed for the random positions and velocity of the particle
 * @param dt Size of timestep in simulation
 * @param totalTime The amount of time the simulation is simulated
 * @param particleInteractions Turns on interactions between particles (True is ON)
**/
void particleSimulationRandom(int seed, double dt, double totalTime, bool particleInteractions){
    double q = 2;
    double m = 40.07;

    double B_0 = 1;
    double V_0 = 10;
    double d = pow(10, 4);

    arma::arma_rng::set_seed(seed);


    //Partikkel 1
    arma::vec r1 = arma::vec(3).randn() * 0.1 * d;  // random initial position
    arma::vec v1 = arma::vec(3).randn() * 0.1 * d;  // random initial velocity


    //Partikkel 2
    arma::vec r2 = arma::vec(3).randn() * 0.1 * d;  // random initial position
    arma::vec v2 = arma::vec(3).randn() * 0.1 * d;  // random initial velocity

    Particle particle_1 = Particle(q, m, r1, v1);
    Particle particle_2 = Particle(q, m, r2, v2);

    std::vector<Particle> particles_in = {particle_1, particle_2};

    PenningTrap trapRK = PenningTrap(B_0, V_0, d, particles_in, particleInteractions);


    double delta = dt;
    int time_steps = totalTime/delta;
    double currentTime = 0;
    int j = 0;

    arma::vec timestamp = arma::vec(time_steps);

    arma::vec posRK_x1 = arma::vec(time_steps);
    arma::vec posRK_y1 = arma::vec(time_steps);
    arma::vec posRK_z1 = arma::vec(time_steps);


    arma::vec posRK_x2 = arma::vec(time_steps);
    arma::vec posRK_y2 = arma::vec(time_steps);
    arma::vec posRK_z2 = arma::vec(time_steps);


    for (int i = 0; i < time_steps; i++){
        j += 1;
        currentTime += delta;

        posRK_x1[i] = trapRK.particles[0].r(0);
        posRK_y1[i] = trapRK.particles[0].r(1);
        posRK_z1[i] = trapRK.particles[0].r(2);

        posRK_x2[i] = trapRK.particles[1].r(0);
        posRK_y2[i] = trapRK.particles[1].r(1);
        posRK_z2[i] = trapRK.particles[1].r(2);


        trapRK.evolve_RK4(delta);

        timestamp[i] = currentTime;
    }

    std::cout << "Finished RK with " << j << " timsteps, total time: "
                    << delta*time_steps << " MicroSeconds..."<<std::endl;

    print_3D_vec(posRK_x1, posRK_y1, posRK_z1, timestamp, time_steps,
                        "positions of a single particle", "posRK1_without.txt", 16);
    print_3D_vec(posRK_x2, posRK_y2, posRK_z2, timestamp, time_steps,
                    "positions of a single particle", "posRK2_without.txt", 16);
}



/**
 * Simulates specific particles in the penning trap
 * @param particles A vector of particles to be put in the penning trap
 * @param dt Size of timestep in simulation
 * @param totalTime The amount of time the simulation is simulated
 * @param particleInteractions Turns on interactions between particles (True is ON)
**/
void particleSimulationSpecific(std::vector<Particle> particles, double dt, double totalTime, bool particleInteractions){
    double B_0 = 1;
    double V_0 = 10;
    double d = pow(10, 4);


    PenningTrap trapRK = PenningTrap(B_0, V_0, d, particles, particleInteractions);


    double delta = dt;
    int time_steps = totalTime/delta;
    double currentTime = 0;
    int j = 0;

    arma::vec timestamp = arma::vec(time_steps);

    arma::vec posRK_x1 = arma::vec(time_steps);
    arma::vec posRK_y1 = arma::vec(time_steps);
    arma::vec posRK_z1 = arma::vec(time_steps);


    arma::vec posRK_x2 = arma::vec(time_steps);
    arma::vec posRK_y2 = arma::vec(time_steps);
    arma::vec posRK_z2 = arma::vec(time_steps);


    for (int i = 0; i < time_steps; i++){
        j += 1;
        currentTime += delta;

        posRK_x1[i] = trapRK.particles[0].r(0);
        posRK_y1[i] = trapRK.particles[0].r(1);
        posRK_z1[i] = trapRK.particles[0].r(2);

        posRK_x2[i] = trapRK.particles[1].r(0);
        posRK_y2[i] = trapRK.particles[1].r(1);
        posRK_z2[i] = trapRK.particles[1].r(2);


        trapRK.evolve_RK4(delta);

        timestamp[i] = currentTime;
    }

    std::cout << "Finished RK with " << j << " timsteps, total time: "
                    << delta*time_steps << " MicroSeconds..."<<std::endl;

    print_3D_vec(posRK_x1, posRK_y1, posRK_z1, timestamp, time_steps,
                        "positions of a single particle", "posRK1_without.txt", 16);
    print_3D_vec(posRK_x2, posRK_y2, posRK_z2, timestamp, time_steps,
                    "positions of a single particle", "posRK2_without.txt", 16);
}
