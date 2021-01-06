/*
Jan,2021 Minju Kum
This program simulates slow-roll eternal inflation with
chaotic potentials, V(phi)=lambda*(phi^n).
Natural units are used; c = hbar = 1
Phi, V(phi) values are normalized with the reduced Plack mass;
 m_Planck = 1/sqrt(8*pi*G)
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
double rng_constant();
double rng_gaussian(double);

int int main(int argc, char const *argv[]) {
  //variables
  long int n_pat, n_step, n;
  double t, dt;
  double phi, dphi_cl, dphi_qu, phi_max, phi_min;
  double V, V_prime;
  double H;
  clock_t start, finish;

  // FILENAME SETTING
  outfile.open("")

  //read in # of runs, # of steps per run, chaotic potential exponent
  cout << "enter # of patches, # of steps, exponent n" << endl;
  cin >> n_pat >> n_step >> n;

  //calculate the constant lambda
  lambda = //A FUNCTION OF n

  start = clock();

  //patch loop
  for (int i = 0; i < n_patch; i++) {
    t = 0;
    phi = rng_constant();       //get phi_initial
    V = lambda*exp(phi,n);
    V_prime = n*lambda*exp(phi,n-1);
    H = sqrt(V/3);
    dt = 1/H;

    writefile(t,phi);

    //evolution for each patch
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
        //HALT AND QUIT
      }

      //calculate V(phi), H(phi)
      V = lambda*exp(phi,n);
      V_prime = n*lambda*exp(phi,n-1);
      H = sqrt(V/3);
      dt = 1/H;

    } //evolution loop

    outfile << endl;

  } //patch loop

  finish = clock();

  return 0;
}

void writefile(double t, double phi)
{
  outfile << setiosflags(ios::showpoint);
  outfile << setw(15) << setprecision(8) << t;
  outfile << setw(15) << setprecision(8) << phi << endl;
}

//random number generator with a uniform real distribution
double rng_constant()
{
  int seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937_64 generator(seed);
  uniform_real_distribution<double> distribution(0.,1.);
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
