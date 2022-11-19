// LEPAGE'S ANALYSIS
// Project for Theoretical and Numerical Aspects of Nuclear Physics
// Bergonzoni Matteo and Tosca Jacopo

// This program makes use of GSL (GNU Scientific Library), for more information and to download it:
// https://www.gnu.org/software/gsl/doc/html/intro.html

// To compile from the terminal use the command:
// g++ phase_shift.cpp -lgsl -lgslcblas -lm

// PHASE SHIFT

#include "definitions.h"

using namespace std ;

// Dirac's delta potential
double dirac (double r, double x_i , double x_f , int discretization ) { 
    double dx = (x_f-x_i)/(DISCRETIZATION-1);
    double height_dirac = 12/dx;
    double distance = r / dx ; // distance in discretization units dx
    return height_dirac * pow(0.5, distance) ; 
    }; 

// short-range potential: a well potential of depth v_depth
double V_s (double r, double V_depth ) {
    if (fabs(r)<1.3) return V_depth;
    else return 0;
    }; 

// smeared delta, to take in account for the UV cutoff, o(a^2)
double delta_a (double r, double cutoff ) {
    return exp(- r*r / (2*cutoff*cutoff) ) / ( pow(2*M_PI,1.5) * pow(cutoff,3) ) ;
    };

// Laplacian of the smeared delta, o(a^4)
double lap_delta (double r, double cutoff ) {
    return exp(- r*r / (2*cutoff*cutoff)) * (r*r / (cutoff * cutoff) - 3) / (pow(2*M_PI,1.5) * pow(cutoff,5));
    };

// effective potential 
// erf() ensures that at r=0 the potential is finite, also for r->0 then V_eff->Coulomb and it introduces a cutoff of order Lambda = 1/a in the momentum space
long double V_eff (double r, double cutoff ) {
    return - alpha * erf( r/(sqrt(2)*cutoff) ) / r + c_1 * pow(cutoff,2) * delta_a(r,a) - d_1 * pow(cutoff,4) * lap_delta(r,cutoff);
    };

// potential 
long double V (double r, double V_depth , double cutoff ) {
    if (!activate_EFT)
        return -alpha/r + V_s(r, V_depth) + activate_dirac * c * dirac(r);
    else
        return V_eff(r, cutoff);
    }; 

// Conversion from energy E to momentum
double k (double x) {
    return (2 * m / (h * h) * (E - V(x)));
    };

// sign function of the product of the two arguments
bool diff_sign (double a, double b) {
    if (a*b >= 0) return false;
    else return true;
    }; 

// mean function
double mean (double a, double b) {
    return (((a+b) / 2.));
    };

// psi_final returns the value of psi_(N-1) given an energy E with left shooting
long double psi_final(double x_i , double x_f , int discretization , double energy , double psi_0 , double s ) {
    Dx = (x_f - x_i)/(discretization -1.);
    psi_j1 = s; // shooting parameter 
    x_j = x_i + Dx;
    for (int i = 1; i < (discretization - 1); i++) {
        // i in [1, discretization-1] because we compute the function from te second point, psi(x_0) e psi(x_1) are already known
        psi_j2 = 2*( (m*Dx*Dx)/(h*h) * (V(x_j) - energy) + 1 )* psi_j1 - psi_j0;
        psi_j0 = psi_j1;
        psi_j1 = psi_j2;
        x_j += Dx; 
        }
    long double fin_val = psi_j2; 
    x_j = 0;
    psi_j0 = 0; 
    psi_j1 = 0; 
    psi_j2 = 0; 
    Dx = 0; 
    return fin_val; 
    };

// psi_mth_point_left computes the value of psi at the m-th step by left shooting
long double psi_mth_point_left (int m_point, double x_i , double x_f , int discretization , double energy , double psi_0 , double s ){
    Dx = (x_f - x_i)/(discretization -1.);  
    psi_j1 = s; // shooting parameter 
    x_j = x_i + Dx;
    E = energy; 
    for (int i = 0; i < (m_point -1); i++) {
        psi_j2 = 2*( (m*Dx*Dx)/(h*h) * (V(x_j) - energy) + 1 )* psi_j1 - psi_j0;
        psi_j0 = psi_j1; psi_j1 = psi_j2; x_j += Dx; 
        }
    long double fin_val = psi_j2; 
    x_j = 0;
    psi_j0 = 0; 
    psi_j1 = 0; 
    psi_j2 = 0; 
    Dx = 0; 
    if (m_point == 1) return s;
    return fin_val; 
    };

// sch_nodes returns the number of nodes counted
int sch_nodes (double x_i, double x_f, int discretization, double energy, double psi_0 , double s ) {
    int nodes_counted = 0;  
    Dx = (x_f - x_i)/(discretization -1.); 
    psi_j1 = s; // shooting parameter 
    E = energy;  
    x_j = x_i + Dx; 
    for (int i = 1; i < (discretization - 1); i++) {
        psi_j2 = 2*((m*Dx*Dx) / (h*h) * (V(x_j) - E) + 1)*psi_j1 - psi_j0;
        if (diff_sign(psi_j1,psi_j2)) nodes_counted ++ ; // nodes_counted counter 
        psi_j0 = psi_j1; psi_j1 = psi_j2; x_j += Dx; 
        }
    x_j = 0;
    psi_j0 = 0; 
    psi_j1 = 0; 
    psi_j2 = 0; 
    Dx = 0; 
    return nodes_counted;  
    };

// energy_finder finds the energy relative to the setted nodes_wanted
double energy_finder (double E_min , double E_max , double x_ii , double x_ff, int discretizationn ) {
    x_initial = x_ii;
    x_final = x_ff; 
    DISCRETIZATION = discretizationn; 
    E1 = E_min;
    E2 = E_max;
    int x = 0; 
    while (true) {
        E_mean = mean(E1, E2); 
        nodes = sch_nodes(x_ii,x_ff,discretizationn, E_mean) ; 
        if(nodes == nodes_wanted) break;
        if (nodes > nodes_wanted) {E2 = E_mean;} else E1 = E_mean ;     
        // x++ ; if (x=100) break ; 
        }
    return E_mean ;
    } ;

// step 5 of the algorithm
bool smear (double e1 , double e2, double x_i , double x_f , int DIS  ) {
    E_mean = mean(E1, E2);
    if (psi_final()* psi_final() > 0) E1 = E_mean;
    else E2 = E_mean; 
    if ( (E2 - E1) < 0.00000001) return true;
    else return false; 
    }

// eigenvalue find the eigenvalue
double eigenvalue (double E_min, double E_max , double x_ii, double x_ff, int discretizationn) {
    DISCRETIZATION = discretizationn; 
    energy_finder (E_min, E_max, x_ii, x_ff, discretizationn);
    while(true) {
        if ( smear() )  break ; 
        energy_finder();
        }
    return E_mean; 
    }

// normalization_single_shooting finds the normalization constant for single shooting from left
void normalization_single_shooting (double x_i , double x_f , int DIS , double energy ) {
    integral_psi = 0; 
    double dx = (x_f - x_i)/(DIS -1.); 
    for (int i = 0; i < DIS; i++) {
        integral_psi += dx * psi_mth_point_left(i, x_i, x_f, DIS, energy, 0);   
        } 
    };

// find_c returns the value of the coupling constant c (see Lepage's article for more details, eq. 5) 
double find_c (double E_lowest , int n ) { 
    return ( E_lowest + 0.5/(n*n) ) * n * n * n * sqrt(M_PI); 
    }; 

// eigenvalues finds the first 15 eigenvalues for some potential and write them on an array 
double* eigenvalues (double E_min, double E_max, double x_i , double x_f , int discretization ){
    nodes_wanted = 0; 
    static double energies [15] ;
    double E_minn = 0, E_maxx =0; 
    // Fill the array energies[15] with the first 15 eigenvalues  
    E_minn = E_min;
    E_maxx = E_max;  
    E_mean = mean(E_minn,E_maxx); 
    for (int i = 0; i< 15; i++) {
        energies [i] = eigenvalue (E_minn, E_maxx, x_i, x_f, discretization); // only left shooting
        nodes_wanted = nodes_wanted + 1;
        //std::cout << "Level: " << i+1 << "S, eigenvalue: " << energies[i] << std::endl; 
        }
    return energies;
    };      
    
// psi_boom_min finds when the eigenfunction explodes (because of the absence of right shooting)
// it returns the last point with 0 derivatives before the divergence
// it does not work for the level 2S
double psi_boom_min (double energy , int discretization , double x_i , double x_f ) {
    double psi = 0; // value of psi at the y-th step
    double psi1 = 0; // value of psi at the (y-1)-th step
    double psi2 = 0; // value of psi at the (y-2)-th step
    double x_stop = 0.;
    double dx = (x_f-x_i)/(discretization+1);
    double boom_trigger = 1.;
    for (int y = 0; y < discretization; y ++){
        psi = psi_mth_point_left(y, x_i, x_f, discretization, energy, 0);
        if ((psi-psi1)*(psi1-psi2)<0) x_stop = (y-1)*dx;
        if (abs(psi) > boom_trigger) break;
        psi2 = psi1;
        psi1 = psi;
        }
    return x_stop;
    }  
       
// psi_boom finds when the eigenfunction explodes (because of the absence of right shooting)
// it returns the last point whith |psi| < stop_trigger (a small quantity close to zero)
double psi_boom (double energy , int discretization , double x_i , double x_f ) {
    double psi = 0;
    double x_stop = 0.;
    double dx = (x_f-x_i)/(discretization+1);
    double boom_trigger = 1.;
    double stop_trigger = 0.00000008;
    for (int y = 0; y < discretization; y ++){
        psi = psi_mth_point_left(y, x_i, x_f, discretization, energy, 0);
        if (abs(psi) < stop_trigger) x_stop = y*dx;
        if (abs(psi) > boom_trigger) break;
        }
    return x_stop;
    }
    
// save_psi creates a file and here saves the information about the normalized eigenfunction
void save_psi (double x_i , double x_f , int discretization , double energy ) {
    E_mean = energy;
    x_initial = x_i; 
    x_final = x_f;
    DISCRETIZATION = discretization; 
    double x_m = x_i;
    double dx = (x_f - x_i) / (discretization-1);
    normalization_single_shooting(x_i, x_f, discretization, energy);
    std::fstream psi_file ; 
    psi_file.open ("eigenfunction_file.txt", std::fstream::out | std::fstream::app) ;
    for (int m = 0; m < discretization; m++){
        psi_file << x_m << ' ' << psi_mth_point_left(m, x_i, x_f, discretization, energy, 0)/abs(integral_psi) <<'\n';
        x_m=x_m+dx;
        }
    std::cout << "\nEigenfunction data (with energy = " << energy << " and in the range [" << x_i << ", " << x_f << "]) have been saved on the file 'eigenfunction_file.text'.\n";
    };

// plot_V writes on a file the data about potential's behaviour in the range [r_min,rmax_] for a plot (with n points)
void plot_V (double r_max , double r_min, int n ) {
    double dr = (r_max - r_min) / (n-1);
    // Coulomb potential
    double V_C [n]; // Array in which we store data about the Coulomb potential
    activate_EFT = false;
    activate_dirac = 0;
    c = 0;
    v_depth = 0;
    for (int i = 0; i < n; i++) {
        V_C [i] = V(r_min + dr*i);
        }
    // short-range interaction + Coulomb potential
    double V_sr [n]; // Array in which we store data about the short-range interaction + Coulomb potential
    v_depth = -1.5;
    for (int i = 0; i < n; i++) {
        V_sr [i] = V(r_min + dr*i);
        }
    // naive approximation (Coulomb + Dirac)
    double V_na [n]; // Array in which we store data about the Dirac + Coulomb potential    
    v_depth = 0;
    c = find_c(lowest_energy, 1);
    c = -5;
    activate_dirac = 1;
    for (int i = 0; i < n; i++) {
        V_na [i] = V(r_min + dr*i);
        }
    // effective theory
     double V_eft [n]; // Array in which we store data about the effective theory
     v_depth = 0;
     a = 1; 
     c_1 = -37.05;
     d_1 = 0.4;
     activate_EFT = true;
     for (int i = 0; i < n; i++) {
         V_eft [i] = V(r_min + dr*i);
         }
    // save on file
    std::fstream V_file ; 
    V_file.open ("potential_file.txt", std::fstream::out | std::fstream::trunc) ;
    for (int i = 0; i < n; i++) {
        V_file << (r_min + dr*i) << ' ' << V_C[i] << ' ' << V_sr[i] << ' ' << V_na[i] << ' ' << V_eft[i] << std::endl;
        }
    std::cout << "\nPotential data in the range [" << r_min << ", " << r_max <<"] have been saved on the file 'potential_file.text'\n";
    }

// the phases (delta) and (delta+2*n*pi), with n integer are equivalent, even if they seem different
// phase_corrector corrects this problem writing the input modulo 2*pi in the range [-pi,pi].
double phase_corrector (double f1) {
    double f2 = fmod(f1, (M_PI));
    return f2;
    }
 
// Phase shift, method 2 
// phase_shift = arg(Gamma(1+ik)), exact phase shift for Coulomb potential
// for more details see Zucchini's notes eq. 5.24.4
double phase_shift2 (double energy ) {
    double k = 1/sqrt(2*abs(energy)); //Sommerfeld parameter with hbar = e_1 = e_2 = m=1
    gsl_sf_result lnr; // here is saved the value of log|Gamma(1+ik)|
    gsl_sf_result arg; // here is saved the value of arg(Gamma(1+ik))
    int status = gsl_sf_lngamma_complex_e(1, k, &lnr, &arg); // Gamma(1+ik)
    return arg.val;
    } 
       
// Phase shift, method 3-bis
// method from Zucchini notes, inverting eq. 5.24.3, it should be valid for our case
double phase_shift3_b (double x , double x_i ,  double x_f , int discretization , double energy ) {
    x_initial = x_i; 
    x_final = x_f;
    DISCRETIZATION = discretization;
    double dx = (x_f - x_i) / (discretization - 1);
    double Dx = dx;
    double x_boom = psi_boom(energy, discretization);
    int m_th = floor(x/dx);
    double En = energy;
    normalization_single_shooting(x_i, x_boom, discretization, En); 
    double chi_mth = 0;
    if ( x < x_boom) {
        chi_mth = psi_mth_point_left(m_th, x_i, x_f, discretization, En, 0)/integral_psi;
        }
    else chi_mth = 0;
    double phase = asin(chi_mth/2) - sqrt(2 * abs(En)) * x + log(2 * sqrt(2*abs(En)) * x) / sqrt(2*abs(En));
    return phase_corrector(phase);
    }

// error_phase is the function to be minimized in order to find c_1 and d_1
// we use the phase shifts of energy 10^(-5) and 10^(-10)
// synthetic data are recived as input: ph_vs_1 and ph_vs_2
double error_phase (double c, double d, double ph_vs_1, double ph_vs_2){
    double x_in = 0;
    double x_fin = 50;
    double xxx = 10; 
    // effective data
    activate_EFT = true;
    a = 1;
    v_depth = 0;
    c_1 = c;
    d_1 = d;
    E_mean = 0.00001;
    double ph_EFT_1 = phase_shift3_b(xxx, x_in, x_fin, 10000);
    E_mean = 0.0000000001;
    double ph_EFT_2 = phase_shift3_b(xxx, x_in, x_fin, 10000);
    double mean_error = mean(abs(ph_vs_1-ph_EFT_1),abs(ph_vs_2-ph_EFT_2));
    return mean_error;
    }
 
// fast_error_phase is a faster and less precise function to be minimized in order to find c_1 and d_1
// we use just the phase shifts at energy 10^(-5)
// synthetic data are recived as input: ph_vs_1
double fast_error_phase (double c, double d, double ph_vs_1) {
    double x_in = 0;
    double x_fin = 50;
    double xxx = 10; 
    // effective data
    activate_EFT = true;
    a = 1;
    v_depth = 0;
    c_1 = c;
    d_1 = d;
    E_mean = 0.00001;
    double ph_EFT_1 = phase_shift3_b(xxx, x_in, x_fin, 10000);
    double error = abs(ph_vs_1-ph_EFT_1);
    return error;
    }
    
// p4mean compute the expectation value of the operator p^4
double p4mean (double x_i , double x_f , int DIS , double energy ) {   
    double integral_p4 = 0; 
    double dx = (x_f - x_i)/(DIS -1.);  
    normalization_single_shooting(x_i, x_f, DISCRETIZATION, energy);
    double psi_values [DIS];
    for (int i = 0; i < DIS; i++) psi_values[i] = psi_mth_point_left(i, x_i, x_f, DIS, energy, 0);  
    double psi_4derivatives_values [DIS];
    for (int i = 0; i < DIS; i++) {
        psi_4derivatives_values[i] = ( psi_values[i] - 4 * psi_values[i+1] + 6 * psi_values[i+2] - 4 * psi_values[i+3] + psi_values[i+4] ) / (dx*dx*dx*dx) ;
        }    
    for (int i = 0; i < DIS-5; i++) {
        integral_p4 += 1./(integral_psi*integral_psi) * dx * ( psi_values[i] * (psi_4derivatives_values[i]) ); 
        } 
    return integral_p4; 
    };

// integral_delta3_initializer compute the expectation value of delta_a for levels 10S, 15S and 20S
void integral_delta3_initializer (double x_i , double x_f , double DIS ) {
    double integral_delta3 = 0;  
    for (int i = 9; i < 21 ; i += 5) {
        nodes_wanted = i;
        E_mean = eigenvalue(-9, 10, 0, 900, 10000);
        integral_delta3 = 0;
        double dx = (range_psi_eff[i] - x_i)/(DIS -1.); 
        double x = 0; 
        int j = 0;
        while (x<= range_psi_eff[i]) {
            integral_delta3 = integral_delta3 + 1./(integral_psi_level_eff[i]*integral_psi_level_eff[i])*dx * delta_a(x,1)*pow(psi_mth_point_left(j, x_i, range_psi_eff[i], DIS, E_means_eff[i], 0),2);
            x += dx;
            j += 1;
            }
        C[i] = integral_delta3;
        //std::cout << "\n C[" << i << "] = " << C[i] << std::endl;
        }
    } ;
 
// C_i_inizializer compute the expectation value of delta_a for levels 1S, 2S, 3S, 4S, 5S and 6S
void C_i_inizializer (double x_i , double x_f , double DIS ) {
    double integral_delta3 = 0;  
    for (int i = 0; i < 6 ; i += 1) {
        nodes_wanted = i;
        E_mean = eigenvalue(-9, 10, 0, 900, 10000);
        integral_delta3 = 0; 
        double dx = (range_psi_eff[i] - x_i)/(DIS -1.); 
        double x = 0; 
        int j = 0;
        while (x<= range_psi_eff[i]) {
            integral_delta3 = integral_delta3 + 1./(integral_psi_level_eff[i]*integral_psi_level_eff[i])*dx * delta_a(x,1)*pow(psi_mth_point_left(j, x_i, range_psi_eff[i], DIS, E_means_eff[i], 0),2);
            x += dx;
            j += 1;
            }
        C[i] = integral_delta3;
        std::cout << "\n C[" << i << "] = " << C[i] << std::endl;
        }
    }

// integral_lap_delta3_initializer compute the expectation value of lap_delta for levels 10S, 15S and 20S
void integral_lap_delta3_initializer (double x_i , double x_f , int DIS ) {
    double integral_lap_delta3 = 0; 
    for (int i = 9; i < 21 ; i += 5) {
        nodes_wanted = i;
        E_mean = eigenvalue(-9, 10, 0, 900, 10000);
        integral_lap_delta3 = 0; 
        double dx = (range_psi_eff[i] - x_i)/(DIS -1.); 
        double x = 0; 
        int j = 0; 
        double lap_values [DIS];
        while (x <= range_psi_eff[i]) {
            lap_values[j] = lap_delta(x,1);
            x = x + dx;
            j += 1;
            }
        x = 0;
        j = 0;  
        while (x<= range_psi_eff[i]) {
            integral_lap_delta3 = integral_lap_delta3 + 1. / (integral_psi_level_eff[i]*integral_psi_level_eff[i]) * dx * lap_values[j] * pow(psi_mth_point_left(j, x_i, range_psi_eff[i], DIS, E_means_eff[i], 0),2);
            x += dx;
            j += 1; 
            } 
        D[i] = integral_lap_delta3;
        std::cout << "\n D[" << i << "] = " << D[i] << std::endl;
        }
    }

// D_i_inizializer compute the expectation value of delta_a for levels 1S, 2S, 3S, 4S, 5S and 6S    
void D_i_inizializer (double x_i , double x_f , int DIS ) {
    double integral_lap_delta3 = 0; 
    for (int i = 0; i < 6 ; i += 1) {
        nodes_wanted = i;
        E_mean = eigenvalue(-9, 10, 0, 900, 10000);
        integral_lap_delta3 = 0; 
        double dx = (range_psi_eff[i] - x_i)/(DIS -1.); 
        double x = 0; 
        int j = 0; 
        double lap_values [DIS];
        while (x <= range_psi_eff[i]) {
            lap_values[j] = lap_delta(x,1);
            x = x + dx;
            j += 1;
            }
        x = 0;
        j = 0;  
        while (x<= range_psi_eff[i]) {
            integral_lap_delta3 = integral_lap_delta3 + 1. / (integral_psi_level_eff[i]*integral_psi_level_eff[i]) * dx * lap_values[j] * pow(psi_mth_point_left(j, x_i, range_psi_eff[i], DIS, E_means_eff[i], 0),2);
            x += dx;
            j += 1; 
            } 
        D[i] = integral_lap_delta3;
        std::cout << "\n D[" << i << "] = " << D[i] << std::endl;
        }
    }    

// system_solver solves the system that fixes the three above parameters
void system_solver () {
    double det_M_inv = 1. / ( p4_eff_10*(C[14]*D[19]-D[14]*C[19]) - C[9]*(p4_eff_15 * D[19] - D[14]*p4_eff_20) + D[9] * (p4_eff_15*C[19] - C[14]*p4_eff_20)) ; 
    Z = det_M_inv * ( (C[14]*D[19]-D[14]*C[19])*p4_true_10 + (D[9]*C[19]-C[9]*D[19])*p4_true_15 + (C[9]*D[14]-D[9]*C[14])*p4_true_20 );
    Gamma = det_M_inv * ( (D[14]*p4_eff_20 - p4_eff_15*D[19])*p4_true_10 + (p4_eff_10*D[19] - D[9]*p4_eff_20)*p4_true_15 + (D[9]*p4_eff_15 - p4_eff_10*D[14])*p4_true_20 );
    eta = det_M_inv * ( (C[19]*p4_eff_15 - p4_eff_20*C[14])*p4_true_10 + (p4_eff_20*C[9] - C[19]*p4_eff_10)*p4_true_15 + (C[14]*p4_eff_10 - p4_eff_15*C[9])*p4_true_20 );
    }

// save_psi_boom saves the range and the normalization of the eigenfunctions(true and effective), and the eigenvalues (true and effective) of the first 20 energy levels
// useful to speed up to not repeat some long calculations
void save_psi_boom () {
    DISCRETIZATION = 10000;
    int N_en_lev = 20;
    // True value of psi_boom (wave function with Coulomb + short-range potential)
    activate_EFT = false;
    activate_dirac = 0;
    c = 0;
    v_depth = -1.5;
    double psi_boom_true[N_en_lev];
    double E_trues[N_en_lev];
    double E_effs[N_en_lev];
    double normalization_trues[N_en_lev];
    double normalization_eff[N_en_lev];
    for (int i = 0; i < N_en_lev; i++) {
        nodes_wanted = i;
        E_mean = eigenvalue(-9,10,0,300,DISCRETIZATION); 
        //if (i == 1) psi_boom_true[i] = psi_boom(E_mean, DISCRETIZATION, 0, 1000);
        //else 
        psi_boom_true[i] = psi_boom_min(E_mean, DISCRETIZATION, 0, 1000);
        E_trues[i] = eigenvalue(-9,10,0,psi_boom_true[i],DISCRETIZATION); 
        normalization_single_shooting(0, psi_boom_true[i], DISCRETIZATION, E_trues[i]);
        normalization_trues[i] = integral_psi;
        std::cout << "Progress:........." << i*50/(N_en_lev-1) << "%........." << std::endl;
        }
    // Effective value of psi_boom (wave function with effective potential)
    v_depth = 0;
    a = 1; 
    c_1 = -37.05; 
    d_1 = 0.4; 
    activate_EFT = true;
    double psi_boom_eff[N_en_lev];
    for (int i = 0; i < N_en_lev; i++) {
        nodes_wanted = i;
        E_mean = eigenvalue(-9,10,0,500,DISCRETIZATION); 
        //if (i == 1) psi_boom_eff[i] = psi_boom(E_mean, DISCRETIZATION, 0, 1000);
        //else
        psi_boom_eff[i] = psi_boom_min(E_mean, DISCRETIZATION, 0, 1000);

        E_effs[i] = eigenvalue(-9,10,0,psi_boom_eff[i],DISCRETIZATION); 
        normalization_single_shooting(0, psi_boom_eff[i], DISCRETIZATION, E_effs[i]);
        normalization_eff[i] = integral_psi;
        std::cout << "Progress:........." << i*50/(N_en_lev-1)+50 << "%........." << std::endl;
        }
    // Save on file
    std::fstream boom_file ; 
    boom_file.open ("psi_boom_file.txt", std::fstream::out | std::fstream::app) ;
    for (int i = 0; i < N_en_lev; i++){
        // in order : psi_boom, mean energy, normalization (first true then effective)
        boom_file << i << ' ' << psi_boom_true[i] << ' ' << psi_boom_eff[i] << ' ' << E_trues[i] << ' ' << E_effs[i] << ' ' << normalization_trues[i] << ' ' << normalization_eff[i] << std::endl;
        }
    std::cout << "\npsi_boom, eigenfunction and eigenvalue data of the first " << N_en_lev << " have been saved on the file 'psi_boom_file.text'.\n";
    }
        
   
int main () {

 // PHASE SHIFTS
 std::cout << "\n PHASE SHIFT\n";
 
 // general parameters
 DISCRETIZATION = 10000;
 double energy = -0.01; 
 int n_en = 9;
 double in_x = 0;
 double fin_x = 50;
 double xx = 10;
 
 /*
 // Validation of the algorithm
 
 std::cout << "\nCoulomb potential\n" << std::endl;
 // Sets global parameters for the Coulomb potential
 activate_EFT = false;
 activate_dirac = 0;
 c = 0;
 v_depth = 0.;
 
 std::cout << "\nMETHOD 2: delta = arg(Gamma(1+ik))  -  phase_shift2()" << std::endl; 
 for (int i = 1; i < n_en ; i++) {   
     E_mean = energy/(pow(10,i));
     std::cout << "energy: " << E_mean << ", phase shift: " << phase_shift2() << std::endl;
     }
 std::cout << "\nMETHOD 3-bis  -  phase_shift3_b()" << std::endl; 
 for (int i = 1; i < n_en ; i++) {   
     E_mean = energy/(pow(10,i));
     std::cout << "energy: " << E_mean << ", phase shift: " << phase_shift3_b() << std::endl;
     }  
 
 
 

 std::cout << "\nCoulomb + short-range potential (synthetic data)\n" << std::endl;
 // Sets global parameters for the Coulomb + short-range potential
 activate_EFT = false;
 activate_dirac = 0;
 c = 0;
 v_depth = -1.5;
 
 /*    
 std::cout << "\nMETHOD 2: delta = arg(Gamma(1+ik))  -  phase_shift2()" << std::endl; 
 for (int i = 1; i < n_en ; i++) {   
     E_mean = energy/(pow(10,i));
     std::cout << "energy: " << E_mean << ", phase shift: " << phase_shift2() << std::endl;
     }   
/*
 E_mean = 0.00001;
 std::cout << "energy: " << E_mean << ", phase shift: " << phase_shift3_b(xx, in_x, fin_x, 10000) << std::endl;
 E_mean = 0.0000000001;
 std::cout << "energy: " << E_mean << ", phase shift: " << phase_shift3_b(xx, in_x, fin_x, 10000) << std::endl;
 // Sets global parameters for the EFT potential
 activate_EFT = true;
 activate_dirac = 0;
 c = 0;
 v_depth = 0;
 a = 1; 
 c_1 = -40; 
 d_1 = 2.5;
 
 std::cout << "\nEFT potential with c_1 = " << c_1 << " and d_1 = " << d_1 << std::endl;
 std::cout << "\nMETHOD 3-bis: phase_shift3_b()" << std::endl; 
 E_mean = 0.00001;
 std::cout << "energy: " << E_mean << ", phase shift: " << phase_shift3_b(xx, in_x, fin_x, 10000) << std::endl;
 E_mean = 0.0000000001;
 std::cout << "energy: " << E_mean << ", phase shift: " << phase_shift3_b(xx, in_x, fin_x, 10000) << std::endl;
 // Another choice for c_1 and d_1
 c_1 = -41; 
 d_1 = 2.5;
 
 std::cout << "\nEFT potential with c_1 = " << c_1 << " and d_1 = " << d_1 << std::endl;
 std::cout << "\nMETHOD 3-bis: phase_shift3_b()" << std::endl; 
 E_mean = 0.00001;
 std::cout << "energy: " << E_mean << ", phase shift: " << phase_shift3_b(xx, in_x, fin_x, 10000) << std::endl;
 E_mean = 0.0000000001;
 std::cout << "energy: " << E_mean << ", phase shift: " << phase_shift3_b(xx, in_x, fin_x, 10000) << std::endl;
 
 // Another choice for c_1 and d_1
 c_1 = -39; 
 d_1 = 2.5;
 
 std::cout << "\nEFT potential with c_1 = " << c_1 << " and d_1 = " << d_1 << std::endl;
 E_mean = 0.00001;
 std::cout << "energy: " << E_mean << ", phase shift: " << phase_shift3_b(xx, in_x, fin_x, 10000) << std::endl;
 E_mean = 0.0000000001;
 std::cout << "energy: " << E_mean << ", phase shift: " << phase_shift3_b(xx, in_x, fin_x, 10000) << std::endl;
 // Another choice for c_1 and d_1
 c_1 = -40; 
 d_1 = 2.0;
 
 std::cout << "\nEFT potential with c_1 = " << c_1 << " and d_1 = " << d_1 << std::endl;
 E_mean = 0.00001;
 std::cout << "energy: " << E_mean << ", phase shift: " << phase_shift3_b(xx, in_x, fin_x, 10000) << std::endl;
 E_mean = 0.0000000001;
 std::cout << "energy: " << E_mean << ", phase shift: " << phase_shift3_b(xx, in_x, fin_x, 10000) << std::endl;
 
 // Another choice for c_1 and d_1
 c_1 = -40; 
 d_1 = 3;
 
 std::cout << "\nEFT potential with c_1 = " << c_1 << " and d_1 = " << d_1 << std::endl;
 E_mean = 0.00001;
 std::cout << "energy: " << E_mean << ", phase shift: " << phase_shift3_b(xx, in_x, fin_x, 10000) << std::endl;
 E_mean = 0.0000000001;
 std::cout << "energy: " << E_mean << ", phase shift: " << phase_shift3_b(xx, in_x, fin_x, 10000) << std::endl; 
  A QUI */ 
 
 
// minimize err_phase

 cout << "Minimization of the error on the phase\n";

 double c_guess = -37.05; //experience says me that these values are the best ones
 double d_guess = 0.4;
 double c_output = 0;
 double d_output = 0;
 double error = 0.;
 double min_error = 10.;
 double x_in = 0;
 double x_fin = 50;
 double xxx = 10; 
 // synthetic data 
 activate_EFT = false;
 activate_dirac = 0;
 c = 0;
 v_depth = -1.5;
 E_mean = 0.00001;
 double ph_vs_1 = phase_shift3_b(xxx, x_in, x_fin, 10000);
 //E_mean = 0.0000000001;
 //double ph_vs_2 = phase_shift3_b(xxx, x_in, x_fin, 10000);
 
 for (double r = -1; r < 1.5; r=r+0.5){
     for (double s = -0.6; s < 0.7; s=s+0.3){
         //error = error_phase(c_guess + r*2, d_guess + s, ph_vs_1, ph_vs_2);
         error = fast_error_phase(c_guess + r, d_guess + s, ph_vs_1);
         std::cout << "for c = " << c_guess + r << ", d = "<< d_guess + s <<" -> error = "<< error<<std::endl;
         if (error < min_error){
             std::cout << "The best until now!" << std::endl;
             min_error = error;
             c_output = c_guess + r;
             d_output = d_guess + s;
             }
         else std::cout << "We have seen better" << std::endl;
         }
     }
 std::cout<<"min err = " << min_error << ", for c_1 = "<< c_output << ", d_1 = " << d_output << std::endl;

 
 
 
 


 }
