import numpy as np
import matplotlib.pyplot as plt

# Read the data into a NumPy array
ny, nx = 21, 21  # These should match the dimensions used in Fortran
psi = np.loadtxt('psi.dat').reshape((ny, nx))
omega = np.loadtxt('omega.dat').reshape((ny, nx))

# Create the contour plot
plt.figure(1)
cp = plt.contourf(psi, 20)
plt.colorbar(cp)
plt.title('Contour Plot of psi')
plt.xlabel('x')
plt.ylabel('y')

plt.figure(2)
cp = plt.contourf(omega, 20)
plt.colorbar(cp)
plt.title('Contour Plot of omega')
plt.xlabel('x')
plt.ylabel('y')

plt.show()
