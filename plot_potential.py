# Import libraries
import numpy  as np
import matplotlib.pyplot as plt

# Discretization
n = 50

# Load data from file
data = np.loadtxt('potential_file.txt', dtype=float)

# Array with r_i
r = [0 for i in range(n)]

for i in range(n) :
   r[i] = data[i, 0]
   
# Array with Coulomb potential
V_C = [0 for i in range(n)]

for i in range(n) :
   V_C[i] = data[i, 1]

# Array with short range interaction + Coulomb potential
V_sr = [0 for i in range(n)]

for i in range(n) :
   V_sr[i] = data[i, 2]
   
# Array with naive approximation
V_na = [0 for i in range(n)]

for i in range(n) :
   V_na[i] = data[i, 3]
   
# Array with EFT
V_eft = [0 for i in range(n)]

for i in range(n) :
   V_eft[i] = data[i, 4]

# Plot
plt.plot(r, V_C, "-b", label = "Coulomb")
plt.plot(r, V_sr, "-y", label = "Short-range + Coulomb")
#plt.plot(r, V_na, "-g", label = "Dirac $\delta$+ Coulomb")
plt.plot(r, V_eft, "-r", label = "EFT")

plt.xlabel("r")
plt.ylabel("V(r)")

plt.legend()

plt.show()
