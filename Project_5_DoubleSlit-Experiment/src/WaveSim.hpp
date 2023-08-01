#ifndef WAVESIM
#define WAVESIM

#include <armadillo>

class WaveSim{
    private:
    public:
    WaveSim(int M_in, double h_in);
    int idx(int i, int j);
    arma::mat GenerateA(int r, arma::vec a);
    arma::mat GenerateB(int r, arma::vec a);
    void Generate_AB(double dt);
    void Update_u();
    void initialize_u(double x_c, double y_c, double sigma_x, double sigma_y, double px, double py);
    void initialize_V(double x_thic, double x_pos, double y_dist, double y_length, double v0);
    void tripple_slit(double x_thic, double x_pos, double y_dist, double y_length, double v0);


    double h;
    public:
    int M; 
    int size;
    arma::sp_cx_mat A;
    arma::sp_cx_mat B;
    arma::cx_colvec u;
    arma::cx_mat V;
};


#endif