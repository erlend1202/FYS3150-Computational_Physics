#include <iostream>

#include "utilities.hpp"
#include "list"


void append_col(arma::cx_vec u, int lenght, float t, char const *name_file){
    
    std::ofstream myfile;
    char const *fileLoc = "../textfiles/";
    std::string fileName = fileLoc;
    fileName += name_file;

    if (t == 0){
        myfile.open(fileName);
    } else{
        myfile.open(fileName, std::ios::out | std::ios::app);
    }

    for (size_t i = 0; i < lenght; i ++){
                myfile << " " << u(i);
        }
    myfile << std::endl;
    myfile.close();
}
