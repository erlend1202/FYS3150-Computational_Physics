#include "utilities.hpp"
#include "list"
/*
This function will write out the contents of a 3D array to file.
----------
Parameters
----------
mat a: vector which we want to write to file, will be divided into three columns.
int length: length of our arrays, denoted by n.
char name_a: name of what variables are in a.
char name_file: name of the txt file we want to write to.
int width: the space which we dedicate for each variable. Purely for formating.
*/
void print_3D_vec(arma::vec x, arma::vec y, arma::vec z, arma::vec timestamp,
                    int lenght, char const *name_a, char const *name_file, int width){
    std::ofstream myfile;

    char const *fileLoc = "../textfiles/";
    std::string fileName = fileLoc;
    fileName += name_file;


    myfile.open(fileName);
            myfile << "#" << " " << name_a << std::endl;
            myfile << "#" << " " << "n = " << lenght << std::endl;
            for (size_t i = 0; i < lenght; i++){
                myfile << "  " << x(i)
                << std::setw(width) << y(i)
                << std::setw(width) << z(i)
                << std::setw(width) << timestamp(i)
                << std::endl;
            }
            myfile.close();
}

void print_f_omega(double f, double omega_V, arma::vec num_escaped_particles, int size_omega, int lenght, std::string name_file){
    std::ofstream myfile;

    char const *fileLoc = "../textfiles/";
    std::string fileName = fileLoc;
    fileName += name_file;

    myfile.open(fileName);
            myfile << "#" << " " << "n = " << lenght << std::endl;
            for (size_t i = 0; i < size_omega; i++){
                myfile << "  " << f
                << std::setw(16) << num_escaped_particles(i)
                << std::endl;
            }
            myfile.close();
}
