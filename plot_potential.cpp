// LEPAGE'S ANALYSIS
// Project for Theoretical and Numerical Aspects of Nuclear Physics
// Bergonzoni Matteo and Tosca Jacopo

// This program makes use of GSL (GNU Scientific Library), for more information and to download it:
// https://www.gnu.org/software/gsl/doc/html/intro.html

// To compile from the terminal use the command:
// g++ plot_potential.cpp -lgsl -lgslcblas -lm

// PLOT POTENTIAL

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
  

   
int main () {

// PLOT POTENTIAL
 eigenvalue(-9, 10, 0, 900, 10000);
 plot_V(3, 0.2);
}
