import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

psi_data = pd.read_csv("outputs/test_1/psi_3-4_25-25.csv", header=None)
psi_data = psi_data.dropna(axis=1)
psi_data = psi_data.values

# Calculate the grid spacing assuming a unit square domain
h = 1.0 / (psi_data.shape[0] - 1)

# Initialize velocity arrays
u = np.zeros_like(psi_data)
v = np.zeros_like(psi_data)

# Calculate u and v velocities
u[-1, 1:-1] = 5.0
for j in range(1, psi_data.shape[1] - 1):
    for i in range(1, psi_data.shape[0] - 1):
        u[i, j] = (psi_data[i+1, j] - psi_data[i-1, j]) / (2 * h)
        v[i, j] = -(psi_data[i, j+1] - psi_data[i, j-1]) / (2 * h)

speed = np.sqrt(u**2 + v**2)

x = np.linspace(0, 1, psi_data.shape[1])
y = np.linspace(0, 1, psi_data.shape[0])

# Quiver plot for velocity field
plt.figure()
quiver_plot = plt.quiver(x, y, u, v, speed, scale=1/h, cmap='jet')
plt.colorbar(quiver_plot, label='Velocity m/s')
plt.title("Velocity field")
plt.xlabel("x")
plt.ylabel("y")
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("./outputs/velocity_field.png", dpi=300, bbox_inches='tight')

# Contour plot for u
plt.figure()
contour_plot_u = plt.contourf(x, y, u, cmap='jet', levels=100)
plt.colorbar(contour_plot_u, label='u Velocity')
contour_lines_u = plt.contour(x, y, u, colors='black', levels=20)
plt.clabel(contour_lines_u, inline=1, fontsize=8)
plt.title("Contour Plot for u")
plt.xlabel("x")
plt.ylabel("y")
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("./outputs/u_contour.png", dpi=300, bbox_inches='tight')

# Contour plot for v
plt.figure()
contour_plot_v = plt.contourf(x, y, v, cmap='jet', levels=100)
plt.colorbar(contour_plot_v, label='v Velocity')
contour_lines_v = plt.contour(x, y, v, colors='black', levels=20)
plt.clabel(contour_lines_v, inline=1, fontsize=8)
plt.title("Contour Plot for v")
plt.xlabel("x")
plt.ylabel("y")
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("./outputs/v_contour.png", dpi=300, bbox_inches='tight')

# Contour plot for psi
plt.figure()
contour_plot_psi = plt.contourf(x, y, psi_data, cmap='RdBu', levels=100)
plt.colorbar(contour_plot_psi, label='Streamfunction')
contour_lines_psi = plt.contour(x, y, psi_data, colors='black', levels=20)
plt.clabel(contour_lines_psi, inline=1, fontsize=8)
plt.title("Streamline of fluid inside the cavity")
plt.xlabel("x")
plt.ylabel("y")
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("./outputs/psi_contour.png", dpi=300, bbox_inches='tight')
