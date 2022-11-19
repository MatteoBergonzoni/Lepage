# Import libraries
import numpy  as np
import matplotlib.pyplot as plt
import math

# Number of discretizations
D = 12

# Number of energy levels for discretization
L = 10

# Load data from file
data = np.loadtxt('error_file.txt', dtype=float)

# Array with the discretizations
discretization = [0 for i in range(D)]

for j in range(D) :
  discretization[j] = data[j, 0]

# Array with the values of energy eigenvalues and their relative errors (absolute value)
true_energies = [0 for i in range(L)]
energies = [0 for i in range(L*D)]
energies_err = [0 for i in range(L*D)]

for l in range(L) :
  t_e = - 0.5 / ((l+1) * (l+1))
  true_energies[l] = round(t_e, 6 - int(math.floor(math.log10(abs(t_e)))) - 1) # round to 6 significant digits
  #true_energies[l] = t_e

for d in range(D) :
  for l in range(1, L+1) :
    energies[l-1+d*L] = data[d, l]
  
for d in range(D):  
  for l in range(L) :
    energies_err[l+d*L] = abs((energies[l+d*L]-true_energies[l])/true_energies[l])

# Mean error for each discretization
mean_err = [0 for i in range(D)]

for d in range(D) :
  mean_err[d] = np.mean(energies_err[d*L:(d+1)*L])

# Error of level 1 energy for each discretization
energy_err_1 = [0 for i in range(D)]

for j in range(D) :
  energy_err_1[j] = energies_err[j*L]

# Error of level 3 energy for each discretization
energy_err_3 = [0 for i in range(D)]

for j in range(D) :
  energy_err_3[j] = energies_err[2+j*L]

# Error of level 5 energy for each discretization
energy_err_5 = [0 for i in range(D)]

for j in range(D) :
  energy_err_5[j] = energies_err[4+j*L]
 
# Creating 2 subplots
fig, ax = plt.subplots(1, 2)

# First plot: Discretization - Mean relative error
ax[0].set_title("Mean relative error (first ten energy levels)")
ax[0].plot(discretization, mean_err)
ax[0].set(xlabel = "Discretization", ylabel = "|Mean relative error|")
ax[0].semilogy()

# Second plot: Discretization - Relative error (for 3 energy levels)
ax[1].set_title("Relative errors for three energy levels")
ax[1].plot(discretization, energy_err_1, "-b", label = "Level 1")
ax[1].plot(discretization, energy_err_3, "-r", label = "Level 3")
ax[1].plot(discretization, energy_err_5, "-g", label = "Level 5")
ax[1].set(xlabel = "Discretization", ylabel = "|Relative error|")
ax[1].legend()
ax[1].semilogy()

plt.show()

