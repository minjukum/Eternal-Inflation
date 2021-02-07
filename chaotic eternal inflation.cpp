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

ofstream outfile1;
ofstream outfile2;
void writefile1(double,double);
void writefile2(int);
double rng_uniform(double,double);
double rng_gaussian(double);

int main(int argc, char const *argv[]) {
  //variables
  long int n_run, n_step, n, n_exit;
  double t, dt;
  double phi, dphi_cl, dphi_qu, phi_max, phi_end;
  double lambda, V, V_prime;
  double fraction;
  double H;
  double Omega;
  clock_t start, finish;

  //read in # of runs, # of timesteps per run, chaotic potential exponent, fraction
  //since V(phi_initial) = (fraction) in Planck units, it should be that (fraction)=<1
  cout << "enter # of runs, # of timesteps per run, exponent n, fraction" << endl;
  cin >> n_run >> n_step >> n >> fraction;

  //filename setting
  char* filename1 = new char[100];
  sprintf(filename1, "n=%ld potential w %ld iterations.txt", n, n_run);
  outfile1.open(filename1);

  char* filename2 = new char[100];
  sprintf(filename2, "efolds n=%ld potential w %ld iterations.txt", n, n_run);
  outfile2.open(filename2);


  //set constants

  /*  lambda is set so that it satisfies the scalar fluctuation amplitude.

      phi_max is set so that V(phi_max) = 1 (sub-Plankian energy condition)
      Then H^2(phi_max) = (1/3), which roughly satisfies H^2/m^2 ~< 1
      (slow-roll condition in quantum regime)

      phi_initial is set so that V(phi_initial) = (fraction) (=< 1).

      phi_end is where the slow-roll assumption breaks down: i.e. epsilon(phi_end)=1  */

  n_exit = 0;
  lambda = pow(2*M_PI,2)*3*pow(n,2)/pow((120*n + n*n/2),1+n/2)*2.1*pow(10,-9);
  phi_max = pow(1./lambda,1./(double)n);
  phi_end = n/sqrt(2);

  //precalculations(only for checking)
  cout << "lambda is " << lambda << endl;
  cout << "phiend is " << phi_end << " and phimax is " << phi_max << endl;
  phi = pow(fraction/lambda,1./(double)n);
  V = lambda*pow(phi,n);
  V_prime = n*lambda*pow(phi,n-1);
  H = sqrt(V/3);
  dt = 1/H;
  Omega = 2*M_PI*M_PI*pow(V_prime,2)/pow(V,3);

  cout << "H is " << H << endl;
  cout << "initial phi is " << phi << endl;
  cout << "dphi_cl is " << -V_prime/V << endl;
  cout << "Omega is " << Omega << endl;

  start = clock();

  //run loop
  for (int i = 0; i < n_run; i++) {
    // cout << i+1 << "th loop is running..." << endl;
    t = 0;
    //phi = rng_uniform(phi_end,phi_max);       //get phi_initial
    //phi = pow(lambda,-1/((double)n+2));       //assumed to be eternal but...
    phi = pow(fraction/lambda,1./(double)n);
    V = lambda*pow(phi,n);
    V_prime = n*lambda*pow(phi,n-1);
    H = sqrt(V/3);
    dt = 1/H;

    writefile1(t,phi);

    //evolution for each run
    for (int j = 0; j < n_step; j++) {
      t += dt;
      //calculate the next phi
      dphi_cl = -V_prime/V;
      dphi_qu = rng_gaussian(H);
      phi = phi + dphi_cl + dphi_qu;
      if (phi > phi_max) {
        phi = phi_max; //barrier
        //cout << "just hit the barrier!" << endl;
      }
      if (phi < 0) {
        phi = 0; //negative values are not allowed
      }

      writefile1(t,phi);

      //check if the reheating condition is met
      if (phi < phi_end) {
        n_exit++;
        // cout << "quit!" << endl;
        Omega = 2*M_PI*M_PI*pow(V_prime,2)/pow(V,3);
        // cout << "Omega is " << Omega << endl;
        // cout << "total e-folding is " << j << endl;
        writefile2(j);
        break;
      }

      //calculate V(phi), H(phi)
      V = lambda*pow(phi,n);
      V_prime = n*lambda*pow(phi,n-1);
      H = sqrt(V/3);
      dt = 1/H;

    } //evolution loop

  } //run loop

  finish = clock();
  outfile1.close();
  outfile2.close();

  //print the results
  cout << "Done! It took " << (float(finish-start)/CLOCKS_PER_SEC) << " seconds." << endl;
  cout << n_exit << " out of " << n_run << " exited inflation." << endl;
  cout << "For more information, see the generated .txt file." << endl;

  return 0;
}

void writefile1(double t, double phi)
{
  outfile1 << setiosflags(ios::showpoint);
  outfile1 << setw(15) << setprecision(8) << t;
  outfile1 << setw(15) << setprecision(8) << phi << endl;
}

void writefile2(int j)
{
  outfile2 << setiosflags(ios::showpoint);
  outfile2 << setw(15) << j << endl;
}

//random number generator with a uniform real distribution
double rng_uniform(double phi_end, double phi_max)
{
  int seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937_64 generator(seed);
  uniform_real_distribution<double> distribution(phi_end,phi_max);
  double phi = distribution(generator);
  return phi;
}

//random number generator with a Gaussian distribution
double rng_gaussian(double H)
{
  int seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937_64 generator(seed);
  normal_distribution<double> distribution(0.,H);
  double phi = distribution(generator);
  return phi;
}
