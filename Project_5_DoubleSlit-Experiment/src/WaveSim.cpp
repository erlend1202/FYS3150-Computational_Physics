#include <iostream>

#include "WaveSim.hpp"


/**
 * Initializer for the class WaveSim.
 * ----------
 * Parameters
 * ----------
 * int M_in: Size for M, which determines the size of our vectors and matrices.
 * int h_in: stepsize in the x and y direction.
**/
WaveSim::WaveSim(int M_in, double h_in){
    M = M_in;
    h = h_in;
    size = (M-2)*(M-2);
    u = arma::cx_colvec(size);
    //!!! remove this later when we know what u is supposed to be !!!
    u.fill(1);
}

/**
 * An index converter, so we convert from a 2D perspective to a 1D.
 * ----------
 * Parameters
 * ----------
 * int i: index value in x-direction
 * int j: index value in y-direction
**/
int WaveSim::idx(int i, int j){
    return i*(M-2) + j;
}

/**
 * Test method so we can see if we generate our A matrix correctly.
 * ----------
 * Parameters
 * ----------
 * int r: Values alongside the diagonal
 * arma::vec a: values for the diagonal. 
 **/
arma::mat WaveSim::GenerateA(int r, arma::vec a){
    int size = (M-2)*(M-2);
    int N = (M-2);
    arma::mat A = arma::mat(size, size);
    A.fill(0);
    int periodic = 0;
    //cases for shortest row of r values
    for (int i = 0; i < size-N; i++){
        A(i,i) = a(i);

        A(i,i+1) = -r;
        A(i+1,i) = -r;
        //cases every M-2 instance
        periodic++;
        if (periodic == M-2){
            A(i,i+1) = 0;
            A(i+1,i) = 0;
            periodic = 0;
        }

        A(i+N,i) = -r;
        A(i,i+N) = -r;
    
    }

    //rest cases for longest row of r values
    for (int i = size-N; i < size-1; i ++){
        A(i,i) = a(i);
        A(i,i+1) = -r;
        A(i+1,i) = -r;

        periodic++;
        if (periodic == M-2){
            A(i,i+1) = 0;
            A(i+1,i) = 0;
            periodic = 0;
        }
    }    

    //last element case
    A(size-1,size-1) = a(size-1);
    return A;
}


/**
 * Test method so we can see if we generate our B matrix correctly.
 * ----------
 * Parameters
 * ----------
 * int r: Values alongside the diagonal
 * arma::vec b: values for the diagonal. 
 **/
arma::mat WaveSim::GenerateB(int r, arma::vec b){
    int size = (M-2)*(M-2);
    int N = (M-2);
    arma::mat B = arma::mat(size, size);
    B.fill(0);
    int periodic = 0;
    //cases for shortest row of r values
    for (int i = 0; i < size-N; i++){
        B(i,i) = b(i);

        B(i,i+1) = -r;
        B(i+1,i) = -r;

        //cases every M-2 instance
        periodic++;
        if (periodic == M-2){
            B(i,i+1) = 0;
            B(i+1,i) = 0;
            periodic = 0;
        }

        B(i+N,i) = -r;
        B(i,i+N) = -r;
    }

    //rest cases for longest row of r values
    for (int i = size-N; i < size-1; i ++){
        B(i,i) = b(i);
        B(i,i+1) = -r;
        B(i+1,i) = -r;

        periodic++;
        if (periodic == M-2){
            B(i,i+1) = 0;
            B(i+1,i) = 0;
            periodic = 0;
        }
    }    

    //last element case
    B(size-1,size-1) = b(size-1);
    return B;
}

/**
 * Method for initializing A and B matrices.
 * ----------
 * Parameters
 * ----------
 * double dt: size of our timesteps. 
 **/
void WaveSim::Generate_AB(double dt){

    int size = (M-2)*(M-2);
    int N = (M-2);
    std::complex<double> ii(0,1);

    A = arma::sp_cx_mat(size,size);

    B = arma::sp_cx_mat(size,size);

    int periodic = 0;

    arma::cx_double r = ii*dt/(2.0*h*h);

    //cases for shortest row of r values
    for (int i = 0; i < size-N; i++){

        A(i,i) = 1.0 + 4.0*r + (ii*dt/2.0)*V(i);
        B(i,i) = 1.0 - 4.0*r - (ii*dt/2.0)*V(i);

        A(i,i+1) = -r;
        A(i+1,i) = -r;
        B(i,i+1) = r;
        B(i+1,i) = r;

        periodic++;
        if (periodic == M-2){
            A(i,i+1) = 0;
            A(i+1,i) = 0;
            B(i,i+1) = 0;
            B(i+1,i) = 0;
            periodic = 0;
        }

        A(i+N,i) = -r;
        A(i,i+N) = -r;
        B(i+N,i) = r;
        B(i,i+N) = r;
    }

    //rest cases for longest row of r values
    for (int i = size-N; i < size-1; i ++){
        A(i,i) = 1.0 + 4.0*r + (ii*dt/2.0)*V(i);
        B(i,i) = 1.0 - 4.0*r - (ii*dt/2.0)*V(i);

        A(i,i+1) = -r;
        A(i+1,i) = -r;
        B(i,i+1) = r;
        B(i+1,i) = r;

        periodic++;
        if (periodic == M-2){
            A(i,i+1) = 0;
            A(i+1,i) = 0;
            B(i,i+1) = 0;
            B(i+1,i) = 0;
            periodic = 0;
        }
        
    }    

    //last element case
    A(size-1,size-1) = 1.0 + 4.0*r + (ii*dt/2.0)*V(size-N);
    B(size-1,size-1) = 1.0 - 4.0*r - (ii*dt/2.0)*V(size-N);
}

/**
 * Method for initializing the vector u.
 * Here we use Gaussian wave packets to calculate the values.
 * ----------
 * Parameters
 * ----------
 * double x_c: Coordinates of the centre of the initial wave packet on the x-axis.
 * double y_c: Coordinates of the centre of the initial wave packet on the y-axis.
 * double sigma_x: Initial widths of the wave packet in the x-axis direction.
 * double sigma_y: Initial widths of the wave packet in the y-axis direction.
 * double px: Wave packet momenta in x direction.
 * double py: Wave packet momenta in y direction.
 **/
void WaveSim::initialize_u(double x_c, double y_c, double sigma_x, double sigma_y, double px, double py){
    double x;
    double y;

    int k;
    std::complex<double> ii(0,1);

    //arma::cx_double sum = 0;
    double sum = 0;
    for (int i = 0; i < (M-2); i++){
        for (int j = 0; j < (M-2); j++){
            x = i*h;
            y = j*h;
            k = idx(i,j);
            if (x == 0 || x == 1 || y == 0 || y == 1){
                u(k) = 0;
            }
            else{
                u(k) = exp(-(x-x_c)*(x-x_c)/(2*sigma_x*sigma_x) - (y-y_c)*(y-y_c)/(2*sigma_y*sigma_y) + ii*px*(x-x_c) + ii*py*(y-y_c));
                sum += (u(k).real() * u(k).real()) + (u(k).imag() * u(k).imag());
            }
        }
    }

    //since sum is a very low number, sqrt and then division will give huge calculation errors. 
    sum = sqrt(sum);
    u = u/sum;
    //to get around the calculation errors, and making sure sum = 1, we repeat the process until this is true.
    // while (fabs(sum - 1) > 0.0001){
    //     sum = 0;
    //     for (int i = 0; i < size; i++){
    //         sum += (u(i).real() * u(i).real()) + (u(i).imag() * u(i).imag());
    //     }
    //     sum = sqrt(sum);
    //     u = u/sum;

    // }

    std::complex<double> p = arma::cdot(u, u);
    u = u/(sqrt(p));
    std::cout << "sum = " << sum << std::endl;


}

/**
 * Method for initializing the matrix V.
 * ----------
 * Parameters
 * ----------
 * double x_thic: Thicness of the wall in the x-direction.
 * double x_pos: Position of the start of the wall in x-direction.
 * double y_dist: Distance between the slits. This is given in y-direction.
 * double y_length: Length of each slit opening in y-direction.
 * double v0: value we give to the barrier. This has to be set as a high value.
 **/
void WaveSim::initialize_V(double x_thic, double x_pos, double y_dist, double y_length, double v0){
    V = arma::cx_mat(M-2, M-2);
    V.fill(0);
    int i_length = std::ceil(x_thic/(h));
    int i_pos = x_pos/(h);
    

    int j_length = std::ceil(y_length/(h));
    int j_dist = std::floor(y_dist/(2*h)); 

    int j_dist_extra = 0;
    if (j_dist*2 < j_length){
        j_dist_extra = 1;
    }
    std::cout << j_dist_extra << std::endl;
    std::cout << i_length  << " " << i_pos << " "<< j_length << " "<< j_dist << " "<< std::endl;
    for (int i = i_pos; i < i_pos + i_length; i++){
        for (int k = 0; k < M-2; k++){
            V(k,i) = v0;
        }
        for (int j = 0; j < j_length; j++){
            V(i_pos + j_dist + j + j_dist_extra, i) = 0;
            V(i_pos - j_dist - j, i) = 0;
        }
    }
    V.save("../textfiles/V.bin");
    arma::vec V = V.as_col();
}

/**
 * Method for initializing the matrix V as a tripple-slit.
 * ----------
 * Parameters
 * ----------
 * double x_thic: Thicness of the wall in the x-direction.
 * double x_pos: Position of the start of the wall in x-direction.
 * double y_dist: Distance between the slits. This is given in y-direction.
 * double y_length: Length of each slit opening in y-direction.
 * double v0: value we give to the barrier. This has to be set as a high value.
 **/
void WaveSim::tripple_slit(double x_thic, double x_pos, double y_dist, double y_length, double v0){
    V = arma::cx_mat(M-2, M-2);
    V.fill(0);
    int i_length = std::ceil(x_thic/(h));
    int i_pos = x_pos/(h);
    

    int j_length = std::ceil(y_length/(h));
    int j_dist = std::floor(y_dist/(2*h)); 

    int j_dist_extra = 0;
    if (j_dist*2 < j_length){
        j_dist_extra = 1;
    }
    std::cout << j_dist_extra << std::endl;
    std::cout << i_length  << " " << i_pos << " "<< j_length << " "<< j_dist << " "<< std::endl;
    for (int i = i_pos; i < i_pos + i_length; i++){
        for (int k = 0; k < M-2; k++){
            V(k,i) = v0;
        }
        for (int j = 0; j < j_length; j++){
            //V(i_pos + j_dist + j + j_dist_extra, i) = 0;
            //V(i_pos - j_dist - j, i) = 0;

            V(i_pos - j + j_length/2, i) = 0;
            V(i_pos - j - j_length*3/2, i) = 0;
            V(i_pos + j + j_length*3/2, i) = 0;
        }
    }
    V.save("../textfiles/V.bin");
    arma::vec V = V.as_col();    
}




/**
 * Method for updating u using spsolve from arma::vec.
 **/
void WaveSim::Update_u(){
    arma::cx_colvec b = B*u;
    u = arma::spsolve(A, b);
}