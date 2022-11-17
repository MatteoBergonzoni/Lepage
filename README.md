# Lepage's Analysis

## Project

Project for the exam of Theoretical and Numerical Aspects of Nuclear Physics

Students: Bergonzoni Matteo and Tosca Jacopo

## General structure

The program is composed of the following parts:

- **definitions.h** contains all the definition of global variables and functions
- **error_discretization.cpp** creates a file 'error_file.text' containing the errors on the eigenvalues computed with different discretizations.
- **error_discretization.py** plots the data containted in 'error_file.text'.
- **plot_psi.cpp** creates a file 'eigenfunction_file.text' containing the data about an eigenfunction in presence of a Coulomb + short-range potential. The specific eigenfunction can be selected modifying the value of the variable *nodes_wanted*. It is also possible to modify the potential acting on the variables *activate_EFT* (to activete the effective potential), *activate_dirac* (to add a delta-like potential in x=0), *v_depth* (the depth of the potential well).
- **plot_psi.py** plots the data contained in 'eigenfunction_file.text'.
- **plot_potential.cpp** creates a file 'potential_file.text' containing the data about some potentials (Coulomb, Coulomb+short-range, Coulomb+Dirac and effective) in a given range (the arguments of the function *plot_V()*).
- **plot_potential.py** plots the data contained in 'potential_file.text'.
- **error_eigenvalue.cpp** prints the data about the first 15 energy levels with different potentials (Coulomb, Colomb+short-range, Coulomb+Dirac, effective of order o(a<sup>2</sup>) and effective of order o(a<sup>4</sup>)). These data are also saved on the file 'eigenvalue_file.text'.
- **error_eigenvalue.py** plots the data contained in 'eigenvalue_file.text' and prints the mean relative errors for each potential. 
