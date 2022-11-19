// LEPAGE'S ANALYSIS
// Project for Theoretical and Numerical Aspects of Nuclear Physics
// Bergonzoni Matteo and Tosca Jacopo

// This program makes use of GSL (GNU Scientific Library), for more information and to download it:
// https://www.gnu.org/software/gsl/doc/html/intro.html

// To compile from the terminal use the command:
// g++ p4.cpp -lgsl -lgslcblas -lm

// P4

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

    
// p4mean compute the expectation value of the operator p^4
double p4mean (double x_i , double x_f , int DIS , double energy ) {   
    double integral_p4 = 0; 
    double dx = (x_f - x_i)/(DIS -1.);  
    normalization_single_shooting(x_i, x_f, DISCRETIZATION, energy);
    double psi_values [DIS];
    for (int i = 0; i < DIS; i++) psi_values[i] = psi_mth_point_left(i, x_i, x_f, DIS, energy, 0);  
    double psi_4derivatives_values [DIS];
    for (int i = 0; i < DIS; i++) {
        psi_4derivatives_values[i] = ( psi_values[i] - 4 * psi_values[i+1] 
                                      + 6 * psi_values[i+2] - 4 * psi_values[i+3] 
                                      + psi_values[i+4] ) / (dx*dx*dx*dx) ;
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
        //std::cout << "\n C[" << i << "] = " << C[i] << std::endl;
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
        //std::cout << "\n D[" << i << "] = " << D[i] << std::endl;
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
        //std::cout << "\n D[" << i << "] = " << D[i] << std::endl;
        }
    }    

// system_solver solves the system that fixes the three above parameters
void system_solver () {
    double det_M_inv = 1. / ( p4_eff_10*(C[14]*D[19]-D[14]*C[19]) - C[9]*(p4_eff_15 * D[19] - D[14]*p4_eff_20) + D[9] * (p4_eff_15*C[19] - C[14]*p4_eff_20)) ; 
    Z = det_M_inv * ( (C[14]*D[19]-D[14]*C[19])*p4_true_10 + (D[9]*C[19]-C[9]*D[19])*p4_true_15 + (C[9]*D[14]-D[9]*C[14])*p4_true_20 );
    Gamma = det_M_inv * ( (D[14]*p4_eff_20 - p4_eff_15*D[19])*p4_true_10 + (p4_eff_10*D[19] - D[9]*p4_eff_20)*p4_true_15 + (D[9]*p4_eff_15 - p4_eff_10*D[14])*p4_true_20 );
    eta = det_M_inv * ( (C[19]*p4_eff_15 - p4_eff_20*C[14])*p4_true_10 + (p4_eff_20*C[9] - C[19]*p4_eff_10)*p4_true_15 + (C[14]*p4_eff_10 - p4_eff_15*C[9])*p4_true_20 );
    }

        
   
int main () {

// INPUT FROM psi_boom_file
 std::ifstream fin("psi_boom_file.txt");
 double reader; 
 for (int i = 0; i< 20; i++) {
     fin >> reader; 
     fin >> range_psi[i];
     fin >> range_psi_eff[i]; 
     fin >> E_means[i];   
     fin >> E_means_eff[i];
     fin >> integral_psi_level[i];
     fin >> integral_psi_level_eff[i]; 
     }
     

 // Fix Z, Gamma and eta using the levels 10S, 15S e 20S
 
 x_initial = 0;
 x_final = 300;
 DISCRETIZATION = 10000;
 // Coulomb + short range interaction (true)
 activate_EFT = false;
 activate_dirac = 0;
 c = 0;
 v_depth = -1.5;
 nodes_wanted = 9; 
 E_mean = E_means[9];
 p4_true_10 = p4mean(0, range_psi[9], 10000, E_mean);
 std::cout << "p4_true_10 = " << p4_true_10 << std::endl;
 nodes_wanted = 14;
 E_mean = E_means[14];
 p4_true_15 = p4mean(0, range_psi[14], 10000, E_mean);
 std::cout << "p4_true_15 = " << p4_true_15 << std::endl;
 nodes_wanted = 19;
 E_mean = E_means[19];
 p4_true_20 = p4mean(0, range_psi[19], 10000, E_mean);
 std::cout << "p4_true_20 = " << p4_true_20 << std::endl;
 // Effective potential 
 v_depth = 0;
 a = 1; 
 c_1 = -37.05; 
 d_1 = 0.4; 
 activate_EFT = true;
 nodes_wanted = 9; 
 E_mean = E_means_eff[9];
 p4_eff_10 = p4mean(0, range_psi[9], 10000, E_mean);
 std::cout << "p4_eff_10 = " << p4_eff_10 << std::endl;
 nodes_wanted = 14;
 E_mean = E_means_eff[14];
 p4_eff_15 = p4mean(0, range_psi[14], 10000, E_mean);
 std::cout << "p4_eff_15 = " << p4_eff_15 << std::endl;
 nodes_wanted = 19;
 E_mean = E_means_eff[19];
 p4_eff_20 = p4mean(0, range_psi[19], 10000, E_mean);
 std::cout << "p4_eff_20 = " << p4_eff_20 << std::endl; 
 integral_delta3_initializer(); 
 integral_lap_delta3_initializer();
 C_i_inizializer();
 D_i_inizializer(); 
 system_solver();
 std::cout << "\nZ = " << Z << std::endl;
 std::cout << "\nGamma = " << Gamma << std::endl;
 std::cout << "\neta = " << eta << '\n' << std::endl;

/*
 // Result of the previous lines 
 Z = 0.350827;
 Gamma = -3557.78;
 eta = -2390.06;
*/


 // Computation of P4
 
  for (int i = 0; i < 6; i++){
     // Coulomb + short range interaction
     activate_EFT = false;
     activate_dirac = 0;
     c = 0;
     v_depth = -1.5;
     nodes_wanted = i; 
     E_mean = eigenvalue(-9,10,0,300,10000) ;
     std::cout << i+1 <<"S (true): " << p4mean(0, range_psi[i], 30000, E_mean) << std::endl ; 
     // Effective potential
     v_depth = 0;
     a = 1; 
     c_1 = -37.05; 
     d_1 = 0.4; 
     activate_EFT = true;
     nodes_wanted = i;
     E_mean = eigenvalue(-9,10,0,300,10000);
     double p4_eff = p4mean(0, range_psi[i], 30000, E_mean);
     std::cout << i+1 <<"S (eff): " << p4_eff << std::endl;
     // Correct effective potential
     std::cout << i+1 <<"S (c_eff): " << (Z*p4_eff + Gamma*C[i] + eta*D[i]) << std::endl;  
     }
 


 }
