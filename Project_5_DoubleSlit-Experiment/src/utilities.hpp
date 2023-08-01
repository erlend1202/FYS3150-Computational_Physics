
#ifndef UTILITIES
#define UTILITIES

#include <fstream>
#include <string>
#include <armadillo>
#include <vector>

//needed to compile on ubuntu
#include <sstream>
#include <iomanip>

void append_col(arma::cx_vec u, int lenght, float t, char const *name_file);


#endif