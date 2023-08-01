#include <armadillo>
#include <iostream>

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
        // std::cout << r << std::endl;
    }

    void setPos(arma::vec pos){
        r = pos;
    }

};
