# Import libraries
import numpy  as np
import matplotlib.pyplot as plt

# Discretization
D = 10000

# Load data from file
data = np.loadtxt('eigenfunction_file.txt', dtype=float)

# Array with x_m
x_m = [0 for i in range(D)]

for m in range(D) :
   x_m[m] = data[m, 0]
   
# Array with psi_m
psi_m = [0 for i in range(D)]

for m in range(D) :
   psi_m[m] = data[m, 1]


# Plot
plt.plot(x_m, psi_m)
plt.xlabel("x")
plt.ylabel("$\psi(x)$")

plt.show()
