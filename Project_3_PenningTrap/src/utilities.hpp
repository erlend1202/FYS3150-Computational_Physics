#ifndef UTILITIES
#define UTILITIES

#include <fstream>
#include <string>
#include <armadillo>

//needed to compile on ubuntu
#include <sstream>
#include <iomanip>

void print_3D_vec(arma::vec x, arma::vec y, arma::vec z, arma::vec timestamp,
                    int lenght, char const *name_a, char const *name_file, int width);
void print_f_omega(double f, double omega_V, arma::vec num_escaped_particles,
                    int size_omega, int lenght, std::string name_file);

#endif
