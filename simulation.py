import matplotlib.pyplot as plt
import numpy as np
import math

from numerical_integrators import numerical_integration_methods
from governing_equations import flat_earth_eom


# vehicle
r_sphere_m = 0.08
m_sphere_kg = 5
J_sphere_kgm2 = 0.4 * m_sphere_kg * r_sphere_m**2

# dectionary
amod = {
    "m_kg" : 1, 
    "Jxz_b_kgm2": 0,
    "Jxx_b_kgm2": J_sphere_kgm2,
    "Jyy_b_kgm2": J_sphere_kgm2,
    "Jzz_b_kgm2": J_sphere_kgm2
}

# initialization
u0_bf_mps = 0
v0_bf_mps = 0
w0_vf_mps = 0
p0_bf_mps = 0
q0_bf_mps = 0
r0_vf_mps = 0
phi0_rad = 90 * math.pi/180
theta0_rad = 0 * math.pi/180
psi0_rad = 0
p10_n_m = 0
p20_n_m = 0
p30_n_m = 0

x0 = np.array([
    u0_bf_mps,
    v0_bf_mps,
    w0_vf_mps,
    p0_bf_mps,
    q0_bf_mps,
    r0_vf_mps,
    phi0_rad,
    theta0_rad,
    psi0_rad,
    p10_n_m,
    p20_n_m,
    p30_n_m])

# num parameters
nx0 = x0.size

t0 = 0.0
tf = 10.0
h_s = 0.005

# numericl integration
t = np.arange(t0, tf + h_s, h_s)
x = np.zeros((nx0, t.size))

x[:, 0] = x0

t, x = numerical_integration_methods.forward_euler(flat_earth_eom.flat_earth_eom, t, x, h_s, amod)

# plot data

fig, axes = plt.subplots(2, 4, figsize=(10, 6))
fig.set_facecolor('black')

# velocity x
axes[0,0].plot(t, x[0,:], color='yellow')
axes[0,0].set_xlabel('Time [s]', color='white')
axes[0,0].set_ylabel('u [m/s]', color='white')
axes[0,0].grid(True)
axes[0,0].set_facecolor('black')
axes[0,0].tick_params(colors='white')

# velocity y
axes[0,1].plot(t, x[1,:], color='yellow')
axes[0,1].set_xlabel('Time [s]', color='white')
axes[0,1].set_ylabel('v [m/s]', color='white')
axes[0,1].grid(True)
axes[0,1].set_facecolor('black')
axes[0,1].tick_params(colors='white')

# velocity z
axes[0,2].plot(t, x[2,:], color='yellow')
axes[0,2].set_xlabel('Time [s]', color='white')
axes[0,2].set_ylabel('w [m/s]', color='white')
axes[0,2].grid(True)
axes[0,2].set_facecolor('black')
axes[0,2].tick_params(colors='white')

# roll rate
axes[0,3].plot(t, x[3,:], color='yellow')
axes[0,3].set_xlabel('Time [s]', color='white')
axes[0,3].set_ylabel('p [rad/s]', color='white')
axes[0,3].grid(True)
axes[0,3].set_facecolor('black')
axes[0,3].tick_params(colors='white')

# roll angle
axes[1,0].plot(t, x[6,:], color='yellow')
axes[1,0].set_xlabel('Time [s]', color='white')
axes[1,0].set_ylabel('phi [rad]', color='white')
axes[1,0].grid(True)
axes[1,0].set_facecolor('black')
axes[1,0].tick_params(colors='white')

# pitch rate
axes[1,1].plot(t, x[4,:], color='yellow')
axes[1,1].set_xlabel('Time [s]', color='white')
axes[1,1].set_ylabel('q [rad/s]', color='white')
axes[1,1].grid(True)
axes[1,1].set_facecolor('black')
axes[1,1].tick_params(colors='white')

# pitch angle
axes[1,2].plot(t, x[7,:], color='yellow')
axes[1,2].set_xlabel('Time [s]', color='white')
axes[1,2].set_ylabel('theta [rad]', color='white')
axes[1,2].grid(True)
axes[1,2].set_facecolor('black')
axes[1,2].tick_params(colors='white')

# yaw angle
axes[1,3].plot(t, x[8,:], color='yellow')
axes[1,3].set_xlabel('Time [s]', color='white')
axes[1,3].set_ylabel('psi [rad]', color='white')
axes[1,3].grid(True)
axes[1,3].set_facecolor('black')
axes[1,3].tick_params(colors='white')

fig.tight_layout()
plt.show()

