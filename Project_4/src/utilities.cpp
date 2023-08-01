#include "utilities.hpp"
#include "omp_rng.hpp"
#include "list"

/**
 * Function for making a random matrix with elements being -1 or 1. 
 * Returns the matrix it made.
 * ----------
 * Parameters
 * ----------
 * int L: size of each matrix dimension.
**/
arma::mat makeMatrixRandom(int L){
    int numbers[] = {-1, 1};
    arma::mat matrix(L,L);
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            matrix(i,j) = omp_rng::get_random_int_mag();
        }    
    }
    return matrix;
}

/**
 * Function for making an ordered matrix witch follows a checkerboard pattern with -1 and 1. 
 * Returns the matrix it made.
 * ----------
 * Parameters
 * ----------
 * int L: size of each matrix dimension.
**/
arma::mat makeMatrixOrdered(int L){
    int numbers[] = {-1,1};
    arma::mat matrix(L,L);
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            //matrix(i,j) = (i + j) % 2 == 0 ? numbers[0] : numbers[1];
            matrix(i,j) = 1;
        }    
    }
    return matrix;
}


/**
 * This function will write out the contents of twwo vectors to a .txt file with two columns.
 * ----------
 * Parameters
 * ----------
 * arma::vec x: vector which we want to write to file, will be written to the first column.
 * arma::vec y: vector which we want to write to file, will be written to the second column.
 * int length: length of our vectors, denoted by n.
 * char name_file: name of the txt file we want to write to.
 * int width: the space which we dedicate for each variable. Purely for formating.
**/
void print_2_vecs(arma::vec x, arma::vec y, int lenght, char const *name_file, int width){
    std::ofstream myfile;
    char const *fileLoc = "../textfiles/";
    std::string fileName = fileLoc;
    fileName += name_file;

    myfile.open(fileName);
            myfile << "#" << " " << "n = " << lenght << std::endl;
            for (size_t i = 0; i < lenght; i++){
                myfile << "  " << x(i)
                << std::setw(width) << y(i)
                << std::endl;
            }
            myfile.close();
}


/**
 * This function will write out the contents of three vectors to a .txt file with three columns.
 * ----------
 * Parameters
 * ----------
 * arma::vec x: vector which we want to write to file, will be written to the first column.
 * arma::vec y: vector which we want to write to file, will be written to the second column.
 * arma::vec z: vector which we want to write to file, will be written to the third column.
 * int length: length of our vectors, denoted by n.
 * char name_file: name of the txt file we want to write to.
 * int width: the space which we dedicate for each variable. Purely for formating.
**/
void print_3_vecs(arma::vec x, arma::vec y, arma::vec z, int lenght, char const *name_file, int width){
    std::ofstream myfile;
    char const *fileLoc = "../textfiles/";
    std::string fileName = fileLoc;
    fileName += name_file;

    myfile.open(fileName);
            myfile << "#" << " " << "n = " << lenght << std::endl;
            for (size_t i = 0; i < lenght; i++){
                myfile << "  " << x(i)
                << std::setw(width) << y(i)
                << std::setw(width) << z(i)
                << std::endl;
            }
            myfile.close();
}

/**
 * This function will write the contents of a list of vectors to a .txt file with n columns.
 * ----------
 * Parameters
 * ----------
 * std::vector<arma::vec> data: List of vectors to be written to file
 * int length: length of our vectors, denoted by n.
 * char name_file: name of the txt file we want to write to.
 * int width: the space which we dedicate for each variable. Purely for formating.
**/
void print_n_vecs(std::vector<arma::vec> data, int lenght, char const *name_file, int width){
    std::ofstream myfile;
    char const *fileLoc = "../textfiles/";
    std::string fileName = fileLoc;
    fileName += name_file;

    myfile.open(fileName);
            myfile << "#" << " " << "n = " << lenght << std::endl;
            for (size_t i = 0; i < lenght; i++){
                for (size_t j = 0; j < data.size(); j++){
                    myfile << std::setw(width) << data[j][i];
                }
                myfile << std::endl;
            }
            myfile.close();
}



/**
 * @brief 
 * 
 * @param start Starting point of distributed points
 * @param end   End point of distributed points
 * @param delta Spacing between points
 * @return std::vector<double> A vector of values with delta spacing
 */
std::vector<double> deltaSpace(double start, double end, double delta){
    std::vector<double> points;
    double currentPoint = start;
    points.push_back(start);
    while (currentPoint < end-delta){
        currentPoint += delta;
        points.push_back(currentPoint);
    }

    return points;
}