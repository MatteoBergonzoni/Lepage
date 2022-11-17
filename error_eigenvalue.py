# Import libraries
import numpy  as np
import matplotlib.pyplot as plt

# Number of energy levels
L = 7



# Load data from file
data = np.loadtxt('eigenvalue_file.txt', dtype=float)

# Array with the binding energies
E = [0 for i in range(L)] # binding energies for Coulomb theory 
E_D = [0 for i in range(L)] # binding energies for Coulomb theory + delta function (1st order PT)
E_Vs = [0 for i in range(L)] # binding energies for Coulomb theory + short-range interaction (true values)
E_eff = [0 for i in range(L)] # binding energies for effective potential o(a^4)
E_eff2 = [0 for i in range(L)] # binding energies for effective potential o(a^2)


# Array with the relative errors of binding energy
err_E = [0 for i in range(L)] # relative errors of binding energy for Coulomb theory
err_E_D = [0 for i in range(L)] # relative errors of binding energy for Coulomb theory + delta function (1st order PT)
err_E_eff = [0 for i in range(L)] # relative errors of binding energy for effective potential o(a^4)
err_E_eff2 = [0 for i in range(L)] # relative errors of binding energy for effective potential o(a^2)

for l in range(L) :
  E[l] = abs(data[l,2])
  E_D[l] = abs(data[l,4])
  E_Vs[l] = abs(data[l,3])
  E_eff[l] = abs(data[l,5])
  E_eff2[l] = abs(data[l,6])
  err_E[l] = abs((data[l,2] - data[l,3]) / data[l,3])
  err_E_D[l] = abs((data[l,4] - data[l,3]) / data[l,3])
  err_E_eff[l] = abs((data[l,5] - data[l,3]) / data[l,3])
  err_E_eff2[l] = abs((data[l,6] - data[l,3]) / data[l,3])

err_E_array = np.array(err_E)
err_E_D_array = np.array(err_E_D)
err_E_eff_array = np.array(err_E_eff)
err_E_eff2_array = np.array(err_E_eff2)

print("mean error for Coulomb potential         = ", err_E_array.mean())
print("mean error for Coulomb + delta potential = ", err_E_D_array.mean())
print("mean error for EFT o(a^4)                = ", err_E_eff_array.mean())
print("mean error for EFT o(a^2)                = ", err_E_eff2_array.mean())


# Plot
plt.plot(E, err_E, ":b", label = "Coulomb theory", marker = ".", markersize = 12)
plt.plot(E_D, err_E_D, ":r", label = "Coulomb + $\delta$ theory", marker = ".", markersize = 12)
plt.plot(E_eff, err_E_eff, ":g", label = "Effective theory o($a^4$)", marker = ".", markersize = 12)
plt.plot(E_eff2, err_E_eff2, ":y", label = "Effective theory o($a^2$)", marker = ".", markersize = 12)
plt.title("Relative errors of binding energy")
plt.xlabel("Energy")
plt.xscale("log")
plt.yscale("log")
plt.ylabel("|Relative error|")
plt.legend()

plt.show()

