/*
Jan,2021 Minju Kum
This program simulates slow-roll eternal inflation with
chaotic potentials, V(phi)=lambda*(phi^n).
Natural units are used; c = hbar = 1
Phi, V(phi) values are normalized with the reduced Plack mass;
 m_Planck = 1/sqrt(8*pi*G)

It traces the evolution of phi, counting the number of runs
that exits the inflation.

*/

#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>
#include <fstream>
#include <time.h>
#include <chrono>

using namespace std;

ofstream outfile;
void writefile(double, double);
double rng_uniform();
double rng_gaussian(double);

int int main(int argc, char const *argv[]) {
  //variables
  long int n_run, n_step, n, n_exit;
  double t, dt;
  double phi, dphi_cl, dphi_qu, phi_max, phi_min;
  double lambda, V, V_prime;
  double H;
  clock_t start, finish;

  //read in # of iterations, # of timesteps per iteration, chaotic potential exponent
  cout << "enter # of iterations, # of timesteps per run, exponent n" << endl;
  cin >> n_run >> n_step >> n;

  //filename setting
  char* filename = new char[100];
  sprintf(filename, "n=%ld potential w %ld iterations.txt", n, n_run);
  outfile.open(filename);

  //set constants
  n_exit = 0;
  phi_max = sqrt(8*M_PI);
  phi_min = n/sqrt(2);
  lambda = pow(2*M_PI,2)*3*(n**2)/pow((120*n + n*n/2),1+n/2)*2.1*pow(10,-9);

  start = clock();

  //run loop
  for (int i = 0; i < n_run; i++) {
    t = 0;
    phi = rng_uniform();       //get phi_initial
    V = lambda*exp(phi,n);
    V_prime = n*lambda*exp(phi,n-1);
    H = sqrt(V/3);
    dt = 1/H;

    writefile(t,phi);

    //evolution for each run
    for (int j = 0; j < n_step; j++) {
      t += dt;
      //calculate the next phi
      dphi_cl = -V_prime/V;
      dphi_qu = rng_gaussian(H);
      phi = phi + dphi_cl + dphi_qu;
      if (phi > phi_max) {
        phi = phi_max; //barrier
      }

      writefile(t,phi);

      //check if the reheating condition is met
      if (phi < phi_min) {
        n_exit++;
        break;
      }

      //calculate V(phi), H(phi)
      V = lambda*exp(phi,n);
      V_prime = n*lambda*exp(phi,n-1);
      H = sqrt(V/3);
      dt = 1/H;

    } //evolution loop

    outfile << endl;

  } //run loop

  finish = clock();
  outfile.close

  cout << "Done! It took " << (float(finish-start)/CLOCKS_PER_SEC) << " seconds." << endl;
  cout << sprintf("%ld out of %ld exited inflation.",n_exit,n_run) << endl;
  cout << "For more information, see the generated .txt file." << endl;

  return 0;
}

void writefile(double t, double phi)
{
  outfile << setiosflags(ios::showpoint);
  outfile << setw(15) << setprecision(8) << t;
  outfile << setw(15) << setprecision(8) << phi << endl;
}

//random number generator with a uniform real distribution
double rng_uniform()
{
  int seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937_64 generator(seed);
  uniform_real_distribution<double> distribution(0.,sqrt(8*M_PI));
  double phi = distribution(generator);
  return phi;
}

//random number generator with a Gaussian distribution
double rng_gaussian(double H)
{
  int seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937_64 generator(seed);
  normal_distribution<> distribution{0.,H};
  double phi = distribution(generator);
  return phi;
}
