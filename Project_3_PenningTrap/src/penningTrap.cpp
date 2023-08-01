#include <armadillo>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

#include "utilities.hpp"

class Particle{
    public:
        double q;       //Charge
        double m;       //Mass
        arma::vec r;    //Position
        arma::vec v;    //Velocity

    Particle(double q_in, double m_in, arma::vec r_in, arma::vec v_in){
        q = q_in;
        m = m_in;
        r = r_in;
        v = v_in;
    }

    void accelerate(arma::vec accel, double dt){
        v += accel*dt;
    }

    void move(double dt){
        r += v*dt;
    }

    void setPos(arma::vec pos){
        r = pos;
    }

};


/**
 * Class for simulating a Penning trap using ECM and RK4.
 **/
class PenningTrap{
    private:
        double k_e = 1.38935333 * pow(10, 5);
        double T = 9.64852558*10;
        double V = 9.64852558*pow(10, 7);
        double V_0_div_d2 = 9.64852558;


    public:
        bool particleInteractions = false;

        double B_0;
        double V_0;
        double d;
        std::vector<Particle> particles;
        arma::vec B;

        double B_in_save;
        double V_in_save;

        double omega_V;
        double f;


    /**
     * Function for initializing the Penning trap.
     *
     * @param B_in The magnetic field strength given in Tesla
     * @param V_in The electric potential between the electrodes given in Volt
     * @param d_in The characteristic dimension of the trap given in micrometer
     * @tparam particles A set of Particle objects
    **/
    PenningTrap(double B_in, double V_in, double d_in,
                    std::vector<Particle> particles_in){
        B_0 = B_in*T;
        V_0 = V_in*V;
        d = d_in;
        particles = particles_in;

        B_in_save = B_in;
        V_in_save = V_in;

        B = {0, 0, B_0};
    }
    /**
     * Function for initializing the Penning trap.
     * Here we have overloaded with extra parameters.
     *
     * @tparam particleInteractions (Default is On) Turns particle interactions on or off
    **/
    PenningTrap(double B_in, double V_in, double d_in,
                    std::vector<Particle> particles_in,
                    bool particleInteractions_in){

        B_0 = B_in*T;
        V_0 = V_in*V;
        d = d_in;
        particles = particles_in;

        B_in_save = B_in;
        V_in_save = V_in;

        particleInteractions = particleInteractions_in;
        B = {0, 0, B_0};
    }

    PenningTrap(double B_in, double V_in, double d_in){
        B_0 = B_in;
        V_0 = V_in;
        d = d_in;

        B_in_save = B_in;
        V_in_save = V_in;

        B = {0, 0, B_0};
    }

    /**
     * Function that adds a particle to the Penning Trap
     *
     * @param p_in A Particle object
    **/
    void add_particle(Particle p_in){
        particles.push_back(p_in);
    }

    /**
     * Function that can turn particle interaction on or off.
     *
     * @param val Boolean expretion, set as true or false
    **/
    void set_particle_interactions(bool val){
        particleInteractions = val;
    }

    /**
     * Function that returns the specified particle
     *
     * @param i Index of the desired particle in the trap
     * @returns A Particle Object
    **/
    Particle return_particle(int i){
        return particles[i];
    }

    /**
     * Function that returns all particles in the Penning Trap
     *
     * @returns A vector of Particle Objects
    **/
    std::vector<Particle> return_particle(){
        return particles;
    }

    /**
     * Evaluates the E field forces on a given particle
     *
     * @param p The particle we want to evaluate
     * @returns A 3D vector describing the force on the particle from the E field
    **/
    arma::vec external_E_field(Particle p){
        double omega_0 = p.q*B_0/p.m;
        arma::vec E = arma::vec(3);

        E(0) = omega_0*p.v(1);
        E(1) = -omega_0*p.v(0);
        E(2) = 0;

        return E;
    }

    /**
     * Evaluates the B field forces on a given particle
     *
     * @param p The particle we want to evaluate
     * @returns A 3D vector describing the force on the particle from the E field
    **/
    arma::vec external_B_field(Particle p){

        double omega_z_2 = 2*p.q*V_0_div_d2/(p.m);
        arma::vec B_thisPos = arma::vec(3);

        B_thisPos(0) = 0.5 * omega_z_2 * p.r(0);
        B_thisPos(1) = 0.5 * omega_z_2 * p.r(1);
        B_thisPos(2) = -omega_z_2 * p.r(2);

        return B_thisPos;
    }



    /**
     * Evaluates the E field forces on a given particle
     * Overloaded to be time dependent
     *
     * @param p The particle we want to evaluate
     * @param t Time step which we want to evaluate the force for
     * @returns A 3D vector describing the force on the particle from the E field
    **/
    arma::vec external_B_field(Particle p, double t){

        double time_V0 = time_dep_V0(t, f, omega_V);


        double omega_z_2 = 2*p.q*time_V0/(d*d*p.m);

        arma::vec B_thisPos = arma::vec(3);
        B_thisPos(0) = 0.5 * omega_z_2 * p.r(0);
        B_thisPos(1) = 0.5 * omega_z_2 * p.r(1);
        B_thisPos(2) = -omega_z_2 * p.r(2);
        return B_thisPos;
    }


    /**
     * Calculates the coulomb force between two particles due to their charge
     *
     * @param i The index of the particle the output is meant for
     * @param j The index of the particle affecting the i particle
     *
     * @returns A 3D vector of the force on particle i caused by particle j,
     *          returns a vector of zeroes if particle interactions are off
     *
    **/
    arma::vec force_particle(int i, int j){
        // Checks if particleInteractions is on
        if (particleInteractions){
            Particle iPart = particles[i];
            Particle jPart = particles[j];

            arma::vec force(3);

            double absDist_3 = pow(sqrt(
                ((iPart.r(0) - jPart.r(0)) * (iPart.r(0) - jPart.r(0))) +
                ((iPart.r(1) - jPart.r(1)) * (iPart.r(1) - jPart.r(1))) +
                ((iPart.r(2) - jPart.r(2)) * (iPart.r(2) - jPart.r(2)))), 3);

            force(0) = k_e * iPart.q*jPart.q * (iPart.r(0) - jPart.r(0))/absDist_3;
            force(1) = k_e * iPart.q*jPart.q * (iPart.r(1) - jPart.r(1))/absDist_3;
            force(2) = k_e * iPart.q*jPart.q * (iPart.r(2) - jPart.r(2))/absDist_3;
            return force;
        }

        // Returns vector of zeroes if particleInteractions is off
        else{
            return arma::vec{0, 0, 0};
        }
    }

    /**
     * Calculates the force from both the E-field and, the B-field on a particle i
     * and sums them.
     *
     * @param i The index of the particle in the Penning trap we want to calculate the forces for
     * @returns A 3D vector representing the force from the E- and B-field
     *          for the given particle i
    **/
    arma::vec externalForces(int i){
        arma::vec Bfield = external_B_field(particles[i]);
        arma::vec Efield = external_E_field(particles[i]);

        return Bfield + Efield;
    }


    /**
     * Calculates the force from the E-field, B-field and other particles
     * on a particle i and sums them.
     *
     * @param i The index of the particle in the Penning trap we want to calculate the forces for
     * @returns A 3D vector representing the force from the E- and B-field and other particles.
     *          for the given particle i
    **/
    arma::vec total_force(int i){
        if (arma::norm(particles[i].r) > d){
            return arma::vec(3).fill(0);
        }
        arma::vec particleForces = {0, 0, 0};

        for (int j = 0; j < particles.size(); j++){
            if (j != i){
                particleForces += force_particle(i, j);
            }
        }


        arma::vec lorentzForce = external_B_field(particles[i]);

        arma::vec eFieldForce = external_E_field(particles[i]);

        // return particleForces/particles[i].m;
        return eFieldForce + lorentzForce + particleForces/particles[i].m;
        //return particleForces/particles[i].m;
    }

    /**
     * Calculates the force from the E-field, B-field and other particles
     * on a particle i and sums them.
     * Overloaded to take time as an additional parameter
     *
     * @param i The index of the particle in the Penning trap we want to calculate the forces for
     * @param t The time instance which we want to calculate the total force for
     * @returns A 3D vector representing the force from the E- and B-field and other particles.
     *          for the given particle i
    **/
    arma::vec total_force(int i, double t){
        if (arma::norm(particles[i].r) > d){
            return arma::vec(3).fill(0);
        }
        arma::vec particleForces = {0, 0, 0};

        if (particleInteractions == 1){

            for (int j = 0; j < particles.size(); j++){
                if (j != i){
                    particleForces += force_particle(i, j);
                }
            }
        }

        arma::vec lorentzForce = external_B_field(particles[i], t);

        arma::vec eFieldForce = external_E_field(particles[i]);

        return eFieldForce + lorentzForce + particleForces/particles[i].m;
        //return particleForces/particles[i].m;
    }

    /**
     * Function to calculate a time dependent V_0
     *
     * @param t Current time
     * @param f amplitude
     * @param omega_V frequency
    **/
    double time_dep_V0(double t, double f, double omega_V){
        double time_dep_V0 = V_0*(1 + f*cos(omega_V*t));
        return time_dep_V0;
    }

    /**
     * Evolves the simulation one timestep forward using the Euler Cromer Method
     *
     * @param dt The size of the timestep given in microseconds
     *
    **/
    void evolve_forward_euler(double dt){

        //Storing forces on all particles
        arma::vec forceOnPart[particles.size()];
        for (int i = 0; i < particles.size(); i++){
            forceOnPart[i] = total_force(i);
        }

        //Updating particles velocity and position
        for (int i = 0; i < particles.size(); i++){

            particles[i].accelerate(forceOnPart[i], dt);
            particles[i].move(dt);
        }
    }


    /**
     * Evolves the simulation one timestep forward using the Runge Kutta Method
     * Overloaded to take in time dependent f and omega_V
     *
     * @param dt The size of the timestep given in microseconds
     * @param currentTime The current time since start of
     *                      simulation given in microseconds
     * @param f Amplitude
     * @param omega_V Frequency
    **/
     // Overload RK4 to take in time dependent f and omega_V
    void evolve_RK4(double dt, double currentTime, double f, double omega_V){

        //Copying current state into a new PenningTrap object
        PenningTrap K1(B_in_save, V_in_save, d, particles);

        //Making K2, K3 and K4 states
        PenningTrap K2(B_in_save, V_in_save, d, particles);
        PenningTrap K3(B_in_save, V_in_save, d, particles);
        PenningTrap K4(B_in_save, V_in_save, d, particles);

        // Set amplitude and frequency
        K1.f = f;
        K1.omega_V = omega_V;
        K2.f = f;
        K2.omega_V = omega_V;
        K3.f = f;
        K3.omega_V = omega_V;
        K4.f = f;
        K4.omega_V = omega_V;

        double t = currentTime;
        #pragma omp parallel
        {
        #pragma omp for
        for (int i = 0; i < particles.size(); i++){
            K2.particles[i].accelerate(K1.total_force(i, t), dt);
        }
        #pragma omp for
        for (int i = 0; i < particles.size(); i++){
            K2.particles[i].move(dt/2);
        }
        #pragma omp for
        for (int i = 0; i < particles.size(); i++){
            K3.particles[i].accelerate(K2.total_force(i, t + dt/2), dt/2);
        }
        #pragma omp for
        for (int i = 0; i < particles.size(); i++){
            K3.particles[i].move(dt/2);
        }
        #pragma omp for
        for (int i = 0; i < particles.size(); i++){
            K4.particles[i].accelerate(K3.total_force(i, t + dt/2), dt);
        }
        #pragma omp for
        for (int i = 0; i < particles.size(); i++){
            K4.particles[i].move(dt/2);
        }
        #pragma omp for
        for (int i = 0; i < particles.size(); i++){
            particles[i].accelerate((double)1/6 * (K1.total_force(i, t + dt) + 2*K2.total_force(i,  t + dt) + 2* K3.total_force(i, t + dt) + K4.total_force(i, t + dt)), dt);
        }
        #pragma omp for
        for (int i = 0; i < particles.size(); i++){
            particles[i].move(dt);
        }
        }

    }


    /**
     * Evolves the simulation one timestep forward using the Runge Kutta Method
     *
     * @param dt The size of the timestep given in microseconds
    **/
    void evolve_RK4(double dt){

        //Copying current state into a new PenningTrap object
        PenningTrap K1(B_in_save, V_in_save, d, particles);


        //Making K2, K3 and K4 states
        PenningTrap K2(B_in_save, V_in_save, d, particles);
        PenningTrap K3(B_in_save, V_in_save, d, particles);
        PenningTrap K4(B_in_save, V_in_save, d, particles);

        //Updating positions and velocties of particles
        for (int i = 0; i < particles.size(); i++){
            K2.particles[i].accelerate(K1.total_force(i), dt/2);
        }
        for (int i = 0; i < particles.size(); i++){
            K2.particles[i].move(dt/2);
        }
        for (int i = 0; i < particles.size(); i++){
            K3.particles[i].accelerate(K2.total_force(i), dt/2);
        }
        for (int i = 0; i < particles.size(); i++){
            K3.particles[i].move(dt/2);
        }
        for (int i = 0; i < particles.size(); i++){
            K4.particles[i].accelerate(K3.total_force(i), dt);
        }
        for (int i = 0; i < particles.size(); i++){
            K4.particles[i].move(dt);
        }

        //Updating positions and velocties of particles
        for (int i = 0; i < particles.size(); i++){
            particles[i].accelerate((double)1/6 * (K1.total_force(i)
                                                    + 2*K2.total_force(i)
                                                    + 2* K3.total_force(i)
                                                    + K4.total_force(i)), dt);
            particles[i].move(dt);

        }
    }

    /**
     * Calculates the number of escaped particles for an instance of
     * the PenningTrap class.
    **/
    int num_escaped_particles(){
        int escaped_particles = 0;
        for (size_t i = 0; i < particles.size(); i ++){
            if (arma::norm(particles[i].r) > d) escaped_particles ++;
        }
        return escaped_particles;
    }
};
