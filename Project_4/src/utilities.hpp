
#ifndef UTILITIES
#define UTILITIES

#include <fstream>
#include <string>
#include <armadillo>
#include <vector>

//needed to compile on ubuntu
#include <sstream>
#include <iomanip>

arma::mat makeMatrixRandom(int L);
arma::mat makeMatrixOrdered(int L);
void print_2_vecs(arma::vec x, arma::vec y, int lenght, char const *name_file, int width);
void print_3_vecs(arma::vec x, arma::vec y, arma::vec z, int lenght, char const *name_file, int width);
void print_n_vecs(std::vector<arma::vec> data, int length, char const *name_file, int width);

std::vector<double> deltaSpace(double start, double end, double delta);

#endif