import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pytictoc import TicToc
import time

# Particles Data
particles = 1000   # Number of particles, should be 1000
starting_height = 1   # Starting height
d = 1000   # Density of detonation material (Kg/m3), assuming all equal
M = 5 * 10**(-2)   # Total mass (Kg)
drag_eff = 1   # Drag efficiency for drag heating
rag_eff = 0.5   # Radiation efficiency for radiation heating
C = 1900   # Heat capacity ((J/K)-Kg)
E = 1   # Detonation energy (J)
Temp = 1500   # Initial temperature (K)
m = (M / particles) * (1 + (1.9 * np.random.random_sample((particles, 1)) - 0.95))   # Mass of particles (Kg)
R = ((3 * m) / (4 * np.pi * d))**(1/3)   # Radius of circular particle (m)
A = np.pi * (R**(2))   # Circular cross sectional area (m2)
SA = 4 * np.pi * (R**(2))   # Surface Area (m2)

# Air and Nature Properties
p_air = 1.225   # Density of air (Kg/m3)
u_air = 1.8 * 10**(-5)   # Viscosity of air
v_surr = np.array([[0.1,0,0]])   # Surrounding velocity
C_air = 1.006 * 10**(3)
g = np.array([[0,0,-9.81]])

# Other parameters
t_step = 10**(-2)   # Seconds per step (seconds)
duration = 1.5   # Duration of simulation (seconds)

#Setting particles starting position
r_particle_start = (2 * np.random.random_sample((particles, 3)) + (-1))/10
r_particle_start[:,2] = r_particle_start[:,2] + 10
r_mass_start = (1/np.sum(m)) * np.sum(m * r_particle_start, axis = 0)

#Setting particles starting velocity
v_particle_start = np.sqrt((2 * E)/ np.sum(m)) * ((r_particle_start - r_mass_start) / np.linalg.norm(r_particle_start - r_mass_start, axis = 1)[:, None])

# Setting initial temperature
t_particle = Temp * np.ones((particles,1), dtype = float)

# Defining drag coefficient function
def drag_coeff(Reynolds_number):
    if Reynolds_number <= 1:
        cd = 24/Reynolds_number
    elif Reynolds_number <= 400:
        cd = 24/(Reynolds_number**0.646)
    elif Reynolds_number <= 3 * 10**(5):
        cd = 0.5
    elif Reynolds_number <= 2 * 10**(6):
        cd = 3.66 * 10**(-4) * Reynolds_number**(0.4275)
    else:
        cd = 0.18
    return cd

fig = plt.figure(figsize = [10,10])
ax = fig.add_subplot(111, projection='3d')
ax.scatter3D(r_particle_start[:,0], r_particle_start[:,1], zs = r_particle_start[:,2], zdir='z', s=20, c='Blue', depthshade=True)
ax.set_xlim(-100, 100)
ax.set_ylim(-100, 100)
ax.set_zlim(0, 100)
ax.view_init(20, 30)
plt.show()
