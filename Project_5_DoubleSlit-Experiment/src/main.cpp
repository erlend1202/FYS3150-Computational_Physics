#include <armadillo>
#include <iostream>

#include "WaveSim.hpp"
#include "utilities.hpp"

// testing functions
void print_AB(WaveSim wave);
void test_A_and_B();
void test_V();
void test_print_matrix(arma::cx_mat mat, int size);
void triple_slit();

// task functions
void task_7();
void task_7_no_barrier();
void single_slit();
void task_8();
void tunneling();

// animation functions
void animate_single();
void animate_double();
void animate_triple();

int main() {
  // test_V();
  // triple_slit();
  // single_slit();
  //   task_7();
  // task_7_no_barrier();
  // task_8();
  // tunneling();
  animate_triple();
}

void triple_slit() {
  double h = 0.005;
  double dt = 2.5 * pow(10, -5);
  int M = (int)1 / h + 1;
  int size = (M - 2) * (M - 2);

  double T = 0.008;
  double x_c = 0.25;
  double y_c = 0.5;
  double sigma_x = 0.05;
  double sigma_y = 0.1;
  double p_x = 200;
  double p_y = 0;
  double v_0 = pow(10, 10);

  WaveSim wave = WaveSim(M, h);

  // Setting all values to zero to get no wall
  double x_thic = 0.02;
  double x_pos = 0.5;
  double y_dist = 0.15;
  double y_length = 0.05;

  wave.tripple_slit(x_thic, x_pos, y_dist, y_length, v_0);
  std::cout << "initialized V" << std::endl;


  wave.initialize_u(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
  std::cout << "initialized u" << std::endl;

  wave.Generate_AB(dt);
  std::cout << "generated AB" << std::endl;

  // wave.initialize_V(x_thic, x_pos, y_dist, y_length, v0);

  int timesteps = int(std::round(T / dt));

  std::cout << timesteps << std::endl;
  std::cout << size << std::endl;

  arma::cx_mat results(size, timesteps);

  for (int i = 0; i < timesteps; i++) {
    wave.Update_u();
    results.col(i) = wave.u;
    std::cout << i << std::endl;
  }

  results.save("../textfiles/results_triple_slit.bin");

  append_col(wave.u, size, 0, "test_u.txt");
}

/**
 * Function that will simulate the experiment.
 * Will also calculate how much the total probability deviates from 1.
 * Using double-slit barrier.
 **/
void task_7() {
  double h = 0.005;
  double dt = 2.5 * pow(10, -5);
  double T = 0.008;
  double x_c = 0.25;
  double sigma_x = 0.05;
  double p_x = 200;
  double y_c = 0.5;
  double sigma_y = 0.1;
  double p_y = 0;
  double v_0 = pow(10, 10);

  int M = (int)1 / h + 1;
  int size = (M - 2) * (M - 2);

  WaveSim wave = WaveSim(M, h);

  double x_thic = 0.02;
  double x_pos = 0.5;
  double y_dist = 0.05;
  double y_length = 0.05;

  wave.initialize_V(x_thic, x_pos, y_dist, y_length, v_0);
  std::cout << "initialized V" << std::endl;

  wave.initialize_u(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
  std::cout << "initialized u" << std::endl;

  wave.Generate_AB(dt);
  std::cout << "generated AB" << std::endl;

  // wave.initialize_V(x_thic, x_pos, y_dist, y_length, v0);

  int timesteps = int(std::round(T / dt));

  std::cout << timesteps << std::endl;
  std::cout << size << std::endl;

  arma::cx_mat results(size, timesteps);

  for (int i = 0; i < timesteps; i++) {
    wave.Update_u();
    results.col(i) = wave.u;
    std::cout << i << std::endl;
  }

  // results.save("../textfiles/results_double_slit.bin");
  results.save("../textfiles/results_double_slit.bin");

  append_col(wave.u, size, 0, "test_u.txt");
}

/**
 * Function for making colormap plots of the evolution over time for our system.
 * It will plot for t=0, t=0.001 and t=0.002.
 * Will also plot the Real and Imaginary parts seperatly aswell.
 **/
void task_8() {
  double h = 0.005;
  double dt = 2.5 * pow(10, -5);
  double T = 0.002;
  double x_c = 0.25;
  double sigma_x = 0.05;
  double p_x = 200;
  double y_c = 0.5;
  double sigma_y = 0.2;
  double p_y = 0;
  double v_0 = pow(10, 10);

  int M = (int)1 / h + 1;
  int size = (M - 2) * (M - 2);

  WaveSim wave = WaveSim(M, h);

  double x_thic = 0.02;
  double x_pos = 0.5;
  double y_dist = 0.05;
  double y_length = 0.05;

  wave.initialize_V(x_thic, x_pos, y_dist, y_length, v_0);
  std::cout << "initialized V" << std::endl;

  wave.initialize_u(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
  std::cout << "initialized u" << std::endl;

  wave.Generate_AB(dt);
  std::cout << "generated AB" << std::endl;

  // wave.initialize_V(x_thic, x_pos, y_dist, y_length, v0);

  int timesteps = int(std::round(T / dt));

  std::cout << timesteps << std::endl;
  std::cout << size << std::endl;

  arma::cx_mat results(size, timesteps);

  for (int i = 0; i < timesteps; i++) {
    wave.Update_u();
    results.col(i) = wave.u;
    std::cout << i << std::endl;
  }

  // results.save("../textfiles/results_slit.bin");
  results.save("../textfiles/results_double_slit_task_8.bin");

  append_col(wave.u, size, 0, "test_u.txt");
}

/**
 * Function for making the data used in the animation of the double-slit
 *experiment. The parameters are the same as in task 8, except T which is
 *prolonged to 15 to make a longer animation.
 **/
void animate_double() {
  double h = 0.005;
  double dt = 2.5 * pow(10, -5);
  double T = 0.015;
  double x_c = 0.25;
  double sigma_x = 0.05;
  double p_x = 200;
  double y_c = 0.5;
  double sigma_y = 0.1;
  double p_y = 0;
  double v_0 = pow(10, 10);

  int M = (int)1 / h + 1;
  int size = (M - 2) * (M - 2);

  WaveSim wave = WaveSim(M, h);

  double x_thic = 0.02;
  double x_pos = 0.5;
  double y_dist = 0.05;
  double y_length = 0.05;

  wave.initialize_V(x_thic, x_pos, y_dist, y_length, v_0);
  std::cout << "initialized V" << std::endl;

  wave.initialize_u(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
  std::cout << "initialized u" << std::endl;

  wave.Generate_AB(dt);
  std::cout << "generated AB" << std::endl;

  int timesteps = int(std::round(T / dt));

  std::cout << timesteps << std::endl;
  std::cout << size << std::endl;

  arma::cx_mat results(size, timesteps);

  for (int i = 0; i < timesteps; i++) {
    wave.Update_u();
    results.col(i) = wave.u;
    std::cout << i << std::endl;
  }

  results.save("../textfiles/results_double_slit_animation.bin");
}

/**
 * Function for making the data used in the animation of a single-slit
 * experiment. The parameters are the same as in task 8, except T which is
 * bprolonged to 15 to make a longer animation.
 **/
void animate_single() {
  double h = 0.005;
  double dt = 2.5 * pow(10, -5);
  int M = (int)1 / h + 1;
  int size = (M - 2) * (M - 2);

  double T = 0.015;
  double x_c = 0.25;
  double y_c = 0.5;
  double sigma_x = 0.05;
  double sigma_y = 0.1;
  double p_x = 200;
  double p_y = 0;
  double v_0 = pow(10, 10);

  WaveSim wave = WaveSim(M, h);

  double x_thic = 0.02;
  double x_pos = 0.5;
  double y_dist = 0;  // Single slit
  double y_length = 0.025;

  wave.initialize_V(x_thic, x_pos, y_dist, y_length, v_0);
  std::cout << "initialized V" << std::endl;

  wave.initialize_u(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
  std::cout << "initialized u" << std::endl;

  wave.Generate_AB(dt);
  std::cout << "generated AB" << std::endl;

  // wave.initialize_V(x_thic, x_pos, y_dist, y_length, v0);

  int timesteps = int(std::round(T / dt));

  std::cout << timesteps << std::endl;
  std::cout << size << std::endl;

  arma::cx_mat results(size, timesteps);

  for (int i = 0; i < timesteps; i++) {
    wave.Update_u();
    results.col(i) = wave.u;
    std::cout << i << std::endl;
  }

  // results.save("../textfiles/results_double_slit.bin");
  results.save("../textfiles/results_single_slit_animation.bin");
}

void animate_triple() {
  double h = 0.005;
  double dt = 2.5 * pow(10, -5);
  int M = (int)1 / h + 1;
  int size = (M - 2) * (M - 2);

  double T = 0.015;
  double x_c = 0.25;
  double y_c = 0.5;
  double sigma_x = 0.05;
  double sigma_y = 0.1;
  double p_x = 200;
  double p_y = 0;
  double v_0 = pow(10, 10);

  WaveSim wave = WaveSim(M, h);

  // Setting all values to zero to get no wall
  double x_thic = 0.02;
  double x_pos = 0.5;
  double y_dist = 0.15;
  double y_length = 0.05;

  wave.initialize_V(x_thic, x_pos, y_dist, y_length, v_0);
  std::cout << "initialized V" << std::endl;

  // Code for triple slit
  arma::vec test_ul = arma::vec(2);
  arma::vec test_lr = arma::vec(2);

  test_ul(0) = 0.48;
  test_ul(1) = 0.475;

  test_lr(0) = 0.52;
  test_lr(1) = 0.525;

  wave.removePot(test_ul, test_lr, pow(10, 10));

  wave.initialize_u(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
  std::cout << "initialized u" << std::endl;

  wave.Generate_AB(dt);
  std::cout << "generated AB" << std::endl;

  // wave.initialize_V(x_thic, x_pos, y_dist, y_length, v0);

  int timesteps = int(std::round(T / dt));

  std::cout << timesteps << std::endl;
  std::cout << size << std::endl;

  arma::cx_mat results(size, timesteps);

  for (int i = 0; i < timesteps; i++) {
    wave.Update_u();
    results.col(i) = wave.u;
    std::cout << i << std::endl;
  }
  results.save("../textfiles/results_triple_slit_animation.bin");
}

/**
 * Function for simulating quantum-tunneling by setting the width of the
 *barrier to 1.
 **/
void tunneling() {
  double h = 0.005;
  double dt = 2.5 * pow(10, -5);
  double T = 0.008;
  double x_c = 0.25;
  double sigma_x = 0.05;
  double p_x = 200;
  double y_c = 0.5;
  double sigma_y = 0.1;
  double p_y = 0;
  double v_0 = pow(10, 10);

  int M = (int)1 / h + 1;
  int size = (M - 2) * (M - 2);

  WaveSim wave = WaveSim(M, h);

  double x_thic = 0.01;
  double x_pos = 0.5;
  double y_dist = 0.;
  double y_length = 0.05;

  wave.initialize_V(x_thic, x_pos, y_dist, y_length, v_0);
  std::cout << "initialized V" << std::endl;

  wave.initialize_u(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
  std::cout << "initialized u" << std::endl;

  wave.Generate_AB(dt);
  std::cout << "generated AB" << std::endl;

  // wave.initialize_V(x_thic, x_pos, y_dist, y_length, v0);

  int timesteps = int(std::round(T / dt));

  std::cout << timesteps << std::endl;
  std::cout << size << std::endl;

  arma::cx_mat results(size, timesteps);

  for (int i = 0; i < timesteps; i++) {
    wave.Update_u();
    results.col(i) = wave.u;
    std::cout << i << std::endl;
  }

  // results.save("../textfiles/results_double_slit.bin");
  results.save("results_double_slit.bin");

  append_col(wave.u, size, 0, "test_u.txt");
}

/**
 * Function that will simulate the experiment.
 * Will also calculate how much the total probability deviates from 1.
 * Using no barrier.
 **/
void task_7_no_barrier() {
  double h = 0.005;
  double dt = 2.5 * pow(10, -5);
  double T = 0.008;
  double x_c = 0.25;
  double sigma_x = 0.05;
  double p_x = 200;
  double y_c = 0.5;
  double sigma_y = 0.05;
  double p_y = 0;
  double v_0 = pow(10, 10);

  int M = (int)1 / h + 1;
  int size = (M - 2) * (M - 2);

  WaveSim wave = WaveSim(M, h);

  // Setting all values to zero to get no wall
  double x_thic = 0;
  double x_pos = 0;
  double y_dist = 0;
  double y_length = 0;

  wave.initialize_V(x_thic, x_pos, y_dist, y_length, v_0);
  std::cout << "initialized V" << std::endl;

  wave.initialize_u(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
  std::cout << "initialized u" << std::endl;

  wave.Generate_AB(dt);
  std::cout << "generated AB" << std::endl;

  // wave.initialize_V(x_thic, x_pos, y_dist, y_length, v0);

  int timesteps = int(std::round(T / dt));

  std::cout << timesteps << std::endl;
  std::cout << size << std::endl;

  arma::cx_mat results(size, timesteps);

  for (int i = 0; i < timesteps; i++) {
    wave.Update_u();
    results.col(i) = wave.u;
    std::cout << i << std::endl;
  }

  results.save("../textfiles/results_no_slit.bin");

  append_col(wave.u, size, 0, "test_u.txt");
}

void single_slit() {
  double h = 0.005;
  double dt = 2.5 * pow(10, -5);
  int M = (int)1 / h + 1;
  int size = (M - 2) * (M - 2);

  double T = 0.008;
  double x_c = 0.25;
  double y_c = 0.5;
  double sigma_x = 0.05;
  double sigma_y = 0.1;
  double p_x = 200;
  double p_y = 0;
  double v_0 = pow(10, 10);

  WaveSim wave = WaveSim(M, h);

  double x_thic = 0.02;
  double x_pos = 0.5;
  double y_dist = 0;  // Single slit
  double y_length = 0.025;

  wave.initialize_V(x_thic, x_pos, y_dist, y_length, v_0);
  std::cout << "initialized V" << std::endl;

  wave.initialize_u(x_c, y_c, sigma_x, sigma_y, p_x, p_y);
  std::cout << "initialized u" << std::endl;

  wave.Generate_AB(dt);
  std::cout << "generated AB" << std::endl;

  // wave.initialize_V(x_thic, x_pos, y_dist, y_length, v0);

  int timesteps = int(std::round(T / dt));

  std::cout << timesteps << std::endl;
  std::cout << size << std::endl;

  arma::cx_mat results(size, timesteps);

  for (int i = 0; i < timesteps; i++) {
    wave.Update_u();
    results.col(i) = wave.u;
    std::cout << i << std::endl;
  }

  results.save("../textfiles/results_single_slit.bin");

  append_col(wave.u, size, 0, "test_u.txt");
}

/**
 * Function for testing that our V implementation is correct.
 **/
void test_V() {
  double h = 0.005;
  double dt = 2.5 * pow(10, -5);
  int M = (int)1 / h + 1;
  int size = (M - 2) * (M - 2);

  double T = 0.008;
  double x_c = 0.25;
  double y_c = 0.5;
  double sigma_x = 0.05;
  double sigma_y = 0.1;
  double p_x = 200;
  double p_y = 0;
  double v_0 = pow(10, 10);

  WaveSim wave = WaveSim(M, h);

  double x_thic = 0.02;
  double x_pos = 0.5;
  double y_dist = 0.05;
  double y_length = 0.05;

  wave.initialize_V(x_thic, x_pos, y_dist, y_length, v_0);
  std::cout << "initialized V" << std::endl;
  // test_print_matrix(wave.V, size);
}

void test_A_and_B() {
  int M = 6;
  int size = (M - 2) * (M - 2);
  double h = 1.0 / size;
  WaveSim wave = WaveSim(M, h);
  arma::vec a = arma::vec((M - 2) * (M - 2));
  arma::vec b = arma::vec((M - 2) * (M - 2));

  a.fill(2);
  b.fill(3);

  int r = 1;

  arma::mat A = wave.GenerateA(r, a);
  arma::mat B = wave.GenerateB(r, b);

  std::cout << A << std::endl;
  std::cout << B << std::endl;
}

void test_print_matrix(arma::cx_mat mat, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (mat(i, j) != (arma::cx_double)0) {
        std::cout << std::setprecision(3) << mat(i, j) << "  ";
      } else {
        std::cout << std::setprecision(1) << mat(i, j) << "      ";
      }
    }
    std::cout << '\n';
    std::cout << '\n';
  }
}

void print_AB(WaveSim wave) {
  int size = (wave.M - 2) * (wave.M - 2);

  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (i == j) {
        std::cout << std::setprecision(3) << wave.A(i, j) << "  ";
      } else {
        std::cout << std::setprecision(2) << wave.A(i, j) << "      ";
      }
    }
    std::cout << '\n';
  }

  std::cout << '\n';

  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (i == j) {
        std::cout << std::setprecision(3) << wave.B(i, j) << "  ";
      } else {
        std::cout << std::setprecision(2) << wave.B(i, j) << "      ";
      }
    }
    std::cout << '\n';
  }
}