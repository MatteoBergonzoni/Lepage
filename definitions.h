#include<iostream>
#include<fstream>
#include<cstdlib>
#include <cmath> 
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_complex.h>

// reduced Planck constant in atomic units
long double h = 1;

// energy
double E = 0;

// small coupling constant which multiplies Dirac's delta 
double c = 0; 

// coupling constants of EFT
double c_1 = 0; // the one that multiplies the smeared delta
double d_1 = 0; // the one that multiplies the Laplacian of the the smeared delta

// mass setted to 1
long double m = 1; 

// fine structure constant in atomic units
long double alpha = 1; 

// discretization of the interval on which we solve the Schroedinger equation
int DISCRETIZATION = 1;

// global variables of x_min e x_max, input of energy_finder() 
double x_final; 
double x_initial; 

// activating variables, used to activate different contributions of the potential
int activate_dirac = 0; 
bool activate_EFT = false;

// Dirac's delta potential
double dirac (double , double =  x_initial, double = x_final, int  = DISCRETIZATION);


// depth of the well potential
double v_depth = 0;

// short-range potential: a well potential of depth v_depth
double V_s (double , double  = v_depth);

// cutoff (in configuration space)
long double a = 0; // inizialized to zero, so we consider all the possible momenta (it has to be modified later)


// smeared delta, to take in account for the UV cutoff, o(a^2)
double delta_a (double , double  = a) ;

// Laplacian of the smeared delta, o(a^4)
double lap_delta (double , double  = a) ;

// effective potential 
// erf() ensures that at r=0 the potential is finite, also for r->0 then V_eff->Coulomb and it introduces a cutoff of order Lambda = 1/a in the momentum space
long double V_eff (double , double  = a) ;


// potential 
long double V (double , double  = v_depth, double  = a) ;

// Conversion from energy E to momentum
double k (double ) ;

// useful global variables used by the function sch_nodes()
double x_j = 0; // position of the j-th step
long double psi_j0 = 0; // value of the eigenfunction at the (j-1)-th step
long double psi_j1 = 0; // value of the eigenfunction at the j-th step j
long double psi_j2 = 0; // value of the eigenfunction at the (j+1)-th step
long double Dx = 0; // interval delta_x 
int nodes = 0; 
double E_mean = 0;
double E1; // E1 takes the temporary value of E_min 
double E2; // E2 takes the temporary value of E_max
int match_point; // used for shooting methood 

// sign function of the product of the two arguments
bool diff_sign (double , double ) ;

// mean function
double mean (double , double ) ;

// psi_final returns the value of psi_(N-1) given an energy E with left shooting
long double psi_final(double  = x_initial, double  = x_final, int  = DISCRETIZATION, double  = E_mean, double  = 0, double  = 0.00001) ;

// psi_mth_point_left computes the value of psi at the m-th step by left shooting
long double psi_mth_point_left (int  = (DISCRETIZATION-1), double  = x_initial, double  = x_final, int  = DISCRETIZATION, double  = E_mean, double  = 0, double  = 0.00001);

// sch_nodes returns the number of nodes counted
int sch_nodes (double , double , int , double , double  = 0, double  = 0.00001) ;

// desired number of nodes, and so defines the energy level
int nodes_wanted = 9 ; 

// energy_finder finds the energy relative to the setted nodes_wanted
double energy_finder (double  = E1, double  = E2, double  = x_initial, double  = x_final, int  = DISCRETIZATION) ;

// step 5 of the algorithm
bool smear (double  = E1, double = E2, double  = x_initial, double  = x_final, int  = DISCRETIZATION ) ;
// eigenvalue find the eigenvalue
double eigenvalue (double , double  , double , double , int ) ;

// normalization
double integral_psi = 0; 

// normalization_single_shooting finds the normalization constant for single shooting from left
void normalization_single_shooting (double  = x_initial, double  = x_final, int  = DISCRETIZATION, double  = E_mean) ;

// lowest energy eigenvalue, useful to define c
double lowest_energy = 0;

// find_c returns the value of the coupling constant c (see Lepage's article for more details, eq. 5) 
double find_c (double  = lowest_energy, int  = 1) ;

// eigenvalues finds the first 15 eigenvalues for some potential and write them on an array 
double* eigenvalues (double , double , double  = x_initial, double  = x_final, int  = DISCRETIZATION);

// psi_boom_min finds when the eigenfunction explodes (because of the absence of right shooting)
// it returns the last point with 0 derivatives before the divergence
// it does not work for the level 2S
double psi_boom_min (double  = E_mean, int  = DISCRETIZATION, double  = x_initial, double  = x_final) ;

// psi_boom finds when the eigenfunction explodes (because of the absence of right shooting)
// it returns the last point whith |psi| < stop_trigger (a small quantity close to zero)
double psi_boom (double  = E_mean, int  = DISCRETIZATION, double  = x_initial, double  = x_final);

// save_psi creates a file and here saves the information about the normalized eigenfunction
void save_psi (double  = x_initial, double  = x_final, int  = DISCRETIZATION, double  = E_mean) ;

// plot_V writes on a file the data about potential's behaviour in the range [r_min,rmax_] for a plot (with n points)
void plot_V (double  = 3, double  = 0.1, int  = 50) ;

// the phases (delta) and (delta+2*n*pi), with n integer are equivalent, even if they seem different
// phase_corrector corrects this problem writing the input modulo 2*pi in the range [-pi,pi].
double phase_corrector (double ) ;

// Phase shift, method 2 
// phase_shift = arg(Gamma(1+ik)), exact phase shift for Coulomb potential
// for more details see Zucchini's notes eq. 5.24.4
double phase_shift2 (double  = E_mean) ;
       
// Phase shift, method 3-bis
// method from Zucchini notes, inverting eq. 5.24.3, it should be valid for our case
double phase_shift3_b (double  = 30, double  = 0,  double  = 50, int  = DISCRETIZATION, double  = E_mean);

// error_phase is the function to be minimized in order to find c_1 and d_1
// we use the phase shifts of energy 10^(-5) and 10^(-10)
// synthetic data are recived as input: ph_vs_1 and ph_vs_2
double error_phase (double , double , double , double );

// fast_error_phase is a faster and less precise function to be minimized in order to find c_1 and d_1
// we use just the phase shifts at energy 10^(-5)
// synthetic data are recived as input: ph_vs_1
double fast_error_phase (double , double , double ) ;

// p4mean compute the expectation value of the operator p^4
double p4mean (double  = x_initial, double  = x_final, int  = DISCRETIZATION, double  = E_mean) ;

// expectation values of the operator p^4 over the states 10S, 15S and 20S (true potential case)
double p4_true_10 = 0;
double p4_true_15 = 0;
double p4_true_20 = 0; 

// expectation values of the operator p^4 over the states 4S, 5S and 6S (true case)
double p4_true_4 = 0;
double p4_true_5 = 0;
double p4_true_6 = 0; 

// expectation values of the operator p^4 over the states 10S, 15S and 20S (effective potential case)
double p4_eff_10 = 0;
double p4_eff_15 = 0;
double p4_eff_20 = 0; 

// expectation values of the operator p^4 over the states 4S, 5S and 6S (effective potential case)
double p4_eff_4 = 0;
double p4_eff_5 = 0;
double p4_eff_6 = 0; 

// some data about the first 20 energy levels in the true and effective case
double range_psi[20]; // true psi_boom data
double range_psi_eff[20]; // effective psi_boom data
double E_means[20]; // true energy
double E_means_eff[20]; // effective energy
double integral_psi_level[20]; // true integral [0, x_boom]
double integral_psi_level_eff[20]; // effective integral [0, x_boom]

// arrays for data about the expectations values of delta_a and its laplacian lap_delta  
double C[20];  // delta_a integral 
double D[20];   // lap_delta integral

// integral_delta3_initializer compute the expectation value of delta_a for levels 10S, 15S and 20S
void integral_delta3_initializer (double  = x_initial, double  = x_final, double  = DISCRETIZATION) ;
 
// C_i_inizializer compute the expectation value of delta_a for levels 1S, 2S, 3S, 4S, 5S and 6S
void C_i_inizializer (double  = x_initial, double  = x_final, double  = DISCRETIZATION) ;

// integral_lap_delta3_initializer compute the expectation value of lap_delta for levels 10S, 15S and 20S
void integral_lap_delta3_initializer (double  = x_initial, double  = x_final, int  = DISCRETIZATION) ;

// D_i_inizializer compute the expectation value of delta_a for levels 1S, 2S, 3S, 4S, 5S and 6S    
void D_i_inizializer (double  = x_initial, double  = x_final, int  = DISCRETIZATION) ;

// parameters used to reproduce the correct expectation value of the operator p^4 in the effective theory  
double Z =0; 
double Gamma = 0; 
double eta = 0; 

// system_solver solves the system that fixes the three above parameters
void system_solver () ;

// save__psi_boom saves the range and the normalization of the eigenfunctions(true and effective), and the eigenvalues (true and effective) of the first 20 energy levels
// useful to speed up to not repeat some long calculations
void save_psi_boom (); // parameters used to reproduce the correct expectation value of the operator p^4 in the effective theory

// Fix Z, Gamma and eta using the levels 10S, 15S e 20S
void fix_constants(); 

//print p^4 mean values in the case of coulomb, EFT and finally EFT corrected (eff)
void p4_print (); 
