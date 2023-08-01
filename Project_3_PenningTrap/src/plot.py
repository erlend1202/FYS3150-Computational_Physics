import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

"""
Makes a plot of the z-axis, a 3D plot, and a 2D plot of the xy-plane
for a simulation of the Penning trap for a single particle.
----------
Parameters
----------
file: a .txt file which holds information of a particles x,y,z 
position for every timestep t.
"""
def singlePlot(file):
    x, y, z, t = np.loadtxt(file, unpack=True, usecols=(0,1,2, 3))

    #Z-plot
    fig = plt.figure(figsize=(8,6))
    plt.plot(t, z, label='z plot')
    plt.xlabel(r'time $\mu s$')
    plt.ylabel(r'Posision z-axis')
    plt.savefig("../figures/z_plot.pdf")
    plt.savefig("../figures/z_plot.pdf")
    plt.show()

    #3D plot
    fig = plt.figure(figsize=(8,6))
    ax = plt.axes(projection='3d')
    ax.plot3D(x, y, z, label='part_1')

    plt.legend()

    plt.show()


    #2D plot XY
    fig = plt.figure(figsize=(8,6))
    plt.plot(x, y)
    plt.show()


"""
Makes a plot of the z-axis, a 3D plot, and a 2D plot of the xy-plane
for a simulation of the Penning trap for two particles.
----------
Parameters
----------
file1: a .txt file which holds information of our first particles x,y,z 
position for every timestep t.
file2: a .txt file which holds information of our second particles x,y,z 
position for every timestep t.
"""
def doublePlot(file1, file2, title):
    x1, y1, z1, t = np.loadtxt(file1, unpack=True, usecols=(0,1,2,3))
    x2, y2, z2 = np.loadtxt(file2, unpack=True, usecols=(0,1,2))

    #Plotting z-movement
    plt.plot(t, z1)
    #plt.plot(t, z2)

    #3D plot
    fig = plt.figure(figsize=(8,6))
    ax = plt.axes(projection='3d')
    ax.plot3D(x1, y1, z1, label=r"particle 1")
    ax.plot3D(x2, y2, z2, label=r"particle 2")
    ax.scatter3D(x1[0], y1[0], z1[0])
    ax.scatter3D(x2[0], y2[0], z2[0])

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_zlabel(r'$z$')

    plt.legend()
    #plt.savefig("../figures/2particle_3d_without.pdf")
    #plt.savefig("../figures/2particle_3d_without.pdf")
    plt.savefig("../figures/2particle_3d.pdf")
    plt.savefig("../figures/2particle_3d.pdf")
    plt.show()


    #2D plot xy
    plt.figure(figsize=(8,6))
    plt.plot(x1, y1, label=r"particle 1")
    plt.plot(x2, y2, label=r"particle 2")
    plt.scatter(x1[0], y1[0], color='green')
    plt.scatter(x2[0], y2[0])

    plt.ylabel(r'$y$')
    plt.xlabel(r'$x$')

    plt.legend()
    #plt.savefig("../figures/2particle_xy_without.pdf")
    #plt.savefig("../figures/2particle_xy_without.pdf")
    plt.savefig("../figures/2particle_xy.pdf")
    plt.savefig("../figures/2particle_xy.pdf")

    plt.show()


"""
Makes a plot of the z-axis, a 3D plot, and a 2D plot of the xy-plane
for a simulation of the Penning trap for three particles.
----------
Parameters
----------
file1: a .txt file which holds information of our first particles x,y,z 
position for every timestep t.
file2: a .txt file which holds information of our second particles x,y,z 
position for every timestep t.
file3: a .txt file which holds information of our third particles x,y,z 
position for every timestep t.
"""
def triplePlot(file1, file2, file3, title):
    x1, y1, z1, t = np.loadtxt(file1, unpack=True, usecols=(0,1,2,3))
    x2, y2, z2 = np.loadtxt(file2, unpack=True, usecols=(0,1,2))
    x3, y3, z3 = np.loadtxt(file3, unpack=True, usecols=(0,1,2))

    #Plotting z-movement
    plt.plot(t, z1)
    plt.plot(t, z2)
    plt.plot(t, z3)

    #3D plot
    fig = plt.figure(figsize=(8,6))
    ax = plt.axes(projection='3d')
    ax.plot3D(x1, y1, z1, label=file1)
    ax.plot3D(x2, y2, z2, label=file2)
    ax.plot3D(x3, y3, z3, label=file3)
    ax.scatter3D(x1[0], y1[0], z1[0])
    ax.scatter3D(x2[0], y2[0], z2[0])
    ax.scatter3D(x3[0], y3[0], z3[0])

    ax.set_title(title)
    plt.legend()
    plt.show()


    #2D plot xy
    fig = plt.figure(figsize=(8,6))
    plt.plot(x1, y1, label=file1)
    plt.plot(x2, y2, label=file2)
    plt.plot(x3, y3, label=file3)
    plt.scatter(x1[0], y1[0], color='green')
    plt.scatter(x2[0], y2[0])
    plt.scatter(x3[0], y3[0])

    plt.title(title)

    plt.legend()
    plt.show()

"""
Calculates the analytical solution in xy-plane.
----------
Parameters
----------
t: an array consisting of all the timesteps we want to calculate the analytical solution for
B_0: The magnetic field strength given in Tesla.
V_0: The electric potential between the electrodes given in Volt
d: The characteristic dimension of the trap given in micrometer 
q: The charge of the particle, given in elementary charge
m: The mass of the particle, given in atomic mass
x_0: Initial position for the particle on the x-axis
v_0: Initial velocity of the particle in y-direction
"""
def generalSol(t, q, m, B_0, V_0, d, x_0, v_0):
    w_0 = q*B_0/m
    w_2z = 2*q*V_0/(m*d*d)

    w_plus = (w_0+np.sqrt(w_0*2 - 2*w_2z))/2
    w_minus = (w_0-np.sqrt(w_0*2 - 2*w_2z))/2

    A_plus = (v_0 + w_minus * x_0)/(w_minus - w_plus)
    A_minus = -(v_0 + w_plus * x_0)/(w_minus - w_plus)

    position = A_plus*np.exp((-1J)*w_plus*t) + A_minus*np.exp((-1J)*w_minus*t)
    return position.real, position.imag

"""
Calculates the analytical solution.
----------
Parameters
----------
t: an array consisting of all the timesteps we want to calculate the analytical solution for
B_0: The magnetic field strength given in Tesla.
V_0: The electric potential between the electrodes given in Volt
d: The characteristic dimension of the trap given in micrometer 
q: The charge of the particle, given in elementary charge
m: The mass of the particle, given in atomic mass
x_0: Initial position for the particle on the x-axis
z_0: Initial position for the particle on the z-axis
v_0: Initial velocity of the particle in y-direction
"""
def analytical3D(t, B_0=1, V_0=10, d=1, q=2, m=40.07, x_0=1, z_0=1.0, v_0=0.2):
    k_e = 1.38935333 * 10**5
    T = 9.64852558 * 10
    V = 9.64852558 * 10**7

    B_0 = B_0*T
    V_0 = V_0*V
    d = d*(10**4)

    V_0_div_d2 = 9.64852558

    w_0 = q*B_0/m
    w_z2 = 2*q*V_0_div_d2/m


    w_m = (w_0 - np.sqrt(w_0**2 - 2*w_z2))/2
    w_p = (w_0 + np.sqrt(w_0**2 - 2*w_z2))/2

    A_p = (v_0 + x_0*w_m)/(w_m - w_p)
    A_m = -(v_0 + x_0*w_p)/(w_m - w_p)

    x = np.zeros(len(t))
    y = np.zeros(len(t))
    z = np.zeros(len(t))

    for i in range(len(t)):
        x[i] = A_p*np.cos(w_p * t[i]) + A_m*np.cos(w_m * t[i])
        y[i] = -(A_p*np.sin(w_p * t[i]) + A_m*np.sin(w_m * t[i]))
        z[i] = z_0*np.cos(np.sqrt(w_z2)*t[i])


    return x, y, z

"""
Plots the analytical solution found using function generalSol.
"""
def plotGenSol():
    t = np.linspace(0, 100, 10000)
    x = np.zeros(10000)
    y = np.zeros(10000)

    for i in range(len(t)):
        x[i], y[i] = generalSol(t[i], 2, 40.07, 1, 10, pow(10, 4), 1, 1)

    fig = plt.figure(figsize=(8,6))
    plt.plot(x, y)
    plt.show()

"""
Makes a plot of the phase space, which is the different dimentions x,y,z
plotted against their respective velocities.
----------
Parameters
----------
filer1: a .txt file which holds information of our first particles x,y,z 
position for every timestep t.
filer2: a .txt file which holds information of our second particles x,y,z 
position for every timestep t.
filev1: a .txt file which holds information of our first particles x,y,z 
velocity for every timestep t.
filev2: a .txt file which holds information of our second particles x,y,z 
velocity for every timestep t.
"""
def velPlots(filev1,filev2, filer1, filer2):
    x1, y1, z1, t = np.loadtxt(filer1, unpack=True, usecols=(0,1,2,3))
    x2, y2, z2 = np.loadtxt(filer2, unpack=True, usecols=(0,1,2))
    vx1, vy1, vz1 = np.loadtxt(filev1, unpack=True, usecols=(0,1,2))
    vx2, vy2, vz2 = np.loadtxt(filev2, unpack=True, usecols=(0,1,2))
    fig = plt.figure(figsize=(8,6))
    plt.plot(x1, vx1, label=r"$x_1$ and $v_{x,1}$")
    plt.plot(x2, vx2, label=r"$x_2$ and $v_{x,2}$")
    plt.xlabel("x-position")
    plt.ylabel(r"$v_x$")
    plt.legend()
    plt.savefig("../figures/vel_x_axis.pdf")
    plt.savefig("../figures/vel_x_axis.pdf")
    plt.show()

    fig = plt.figure(figsize=(8,6))
    plt.plot(y1, vy1, label=r"$y_1$ and $v_{y,1}$")
    plt.plot(y2, vy2, label=r"$y_2$ and $v_{y,2}$")
    plt.xlabel("y-position")
    plt.ylabel(r"$v_y$")
    plt.legend()
    plt.savefig("../figures/vel_y_axis.pdf")
    plt.savefig("../figures/vel_y_axis.pdf")
    plt.show()

    fig = plt.figure(figsize=(8,6))
    plt.plot(z1, vz1, label=r"$z_1$ and $v_{z,1}$")
    plt.plot(z2, vz2, label=r"$z_2$ and $v_{z,2}$")
    plt.xlabel("z-position")
    plt.ylabel(r"$v_z$")
    plt.legend()
    plt.savefig("../figures/vel_z_axis.pdf")
    plt.savefig("../figures/vel_z_axis.pdf")
    plt.show()

"""
Makes a subplot of the phase space, which is the different dimentions x,y,z
plotted against their respective velocities. We subplot with and without 
particle interactions against each other.
----------
Parameters
----------
Files with particle interactions.
filer1: a .txt file which holds information of our first particles x,y,z 
position for every timestep t. 
filer2: a .txt file which holds information of our second particles x,y,z 
position for every timestep t.
filev1: a .txt file which holds information of our first particles x,y,z 
velocity for every timestep t.
filev2: a .txt file which holds information of our second particles x,y,z 
velocity for every timestep t.

Files without particle interactions.
filer3: a .txt file which holds information of our first particles x,y,z 
position for every timestep t.
filer4: a .txt file which holds information of our second particles x,y,z 
position for every timestep t.
filev3: a .txt file which holds information of our first particles x,y,z 
velocity for every timestep t.
filev4: a .txt file which holds information of our second particles x,y,z 
velocity for every timestep t.
"""
def vel_SUBPlots(filev1,filev2, filer1, filer2, filev3,filev4, filer3, filer4):
    x1, y1, z1, t = np.loadtxt(filer1, unpack=True, usecols=(0,1,2,3))
    x2, y2, z2 = np.loadtxt(filer2, unpack=True, usecols=(0,1,2))
    vx1, vy1, vz1 = np.loadtxt(filev1, unpack=True, usecols=(0,1,2))
    vx2, vy2, vz2 = np.loadtxt(filev2, unpack=True, usecols=(0,1,2))

    x3, y3, z3, t1 = np.loadtxt(filer3, unpack=True, usecols=(0,1,2,3))
    x4, y4, z4 = np.loadtxt(filer4, unpack=True, usecols=(0,1,2))
    vx3, vy3, vz3 = np.loadtxt(filev3, unpack=True, usecols=(0,1,2))
    vx4, vy4, vz4 = np.loadtxt(filev4, unpack=True, usecols=(0,1,2))


    fig, axs = plt.subplots(2,figsize=(8,6))
    #fig.suptitle("phase space plots of x and vx")

    axs[0].plot(x1, vx1, label=r"$x_1$ and $v_{x_1}$")
    axs[0].plot(x2, vx2, label=r"$x_2$ and $v_{x_2}$")
    axs[1].plot(x3, vx3, label=r"$x_1$ and $v_{x_1}$")
    axs[1].plot(x4, vx4, label=r"$x_2$ and $v_{x_2}$")

    #axs[0].set_title("with collision")
    #axs[1].set_title("without collision")

    fig.text(0.5, 0.04, r'x-position', ha='center')
    fig.text(0.04, 0.5, r"$v_x$", va='center', rotation='vertical')

    axs[0].legend(loc=4)
    axs[1].legend(loc=4)
    plt.savefig("../figures/vel_x_axis.pdf")
    plt.savefig("../figures/vel_x_axis.pdf")
    plt.show()

    fig, axs = plt.subplots(2,figsize=(8,6))
    #fig.suptitle("phase space plots of x and vx")
    axs[0].plot(y1, vy1, label=r"$y_1$ and $v_{y_1}$")
    axs[0].plot(y2, vy2, label=r"$y_2$ and $v_{y_2}$")
    axs[1].plot(y3, vy3, label=r"$y_1$ and $v_{y_1}$")
    axs[1].plot(y4, vy4, label=r"$y_2$ and $v_{y_2}$")

    """axs[0].set_title("with collision")
    axs[1].set_title("without collision")"""

    fig.text(0.5, 0.04, r'y-position', ha='center')
    fig.text(0.04, 0.5, r"$v_y$", va='center', rotation='vertical')

    axs[0].legend(loc=4)
    axs[1].legend(loc=4)
    plt.savefig("../figures/vel_y_axis.pdf")
    plt.savefig("../figures/vel_y_axis.pdf")
    plt.show()

    fig, axs = plt.subplots(2,figsize=(8,6))
    #fig.suptitle("phase space plots of x and vx")
    axs[0].plot(z1, vz1, label=r"$z_1$ and $v_{z_1}$")
    axs[0].plot(z2, vz2, label=r"$z_2$ and $v_{z_2}$")
    axs[1].plot(z3, vz3, label=r"$z_1$ and $v_{z_1}$")
    axs[1].plot(z4, vz4, label=r"$z_2$ and $v_{z_2}$")

    """axs[0].set_title("with collision")
    axs[1].set_title("without collision")"""

    fig.text(0.5, 0.04, r'z-position', ha='center')
    fig.text(0.04, 0.5, r"$v_z$", va='center', rotation='vertical')

    axs[0].legend(loc=4)
    axs[1].legend(loc=4)
    plt.savefig("../figures/vel_z_axis.pdf")
    plt.savefig("../figures/vel_z_axis.pdf")
    plt.show()


"""
Calculates the error between our numerical solution and an 
analytical solution with same initial conditions. 
Returns an array with errors for given time, and the corresponding time 
array.
----------
Parameters
----------
file: a .txt file which holds information of our particles x,y,z 
position for every timestep t. 

"""
def errorArray(file):
    x, y, z, t = np.loadtxt(file, unpack=True, usecols=(0,1,2,3))

    B_0, V_0, d = np.loadtxt("penningTrapConfig.txt", usecols=(0,1,2))
    # B_0 = 1; V_0 = 10; d = 1
    q = 2
    m = 40
    x_0 = 0
    z_0 = 1
    v_0 = 5


    exX, exY, exZ = analytical3D(t, B_0, V_0, d, q, m, x_0, z_0, v_0)

    error = np.copy(t)
    for i in range(len(t)):
        err_x = np.abs(x[i] - exX[i])
        err_y = np.abs(y[i] - exY[i])
        err_z = np.abs(z[i] - exZ[i])
        #dist_est = np.sqrt(x[i]**2 + y[i]**2 + z[i]**2)
        error[i] = np.abs(1-np.sqrt(err_x**2 + err_y**2 + err_z**2))

    return error, t

"""
Plots the error between our numerical solution and an analytical 
solution with same initial conditions for multiple simulations
with different timesteps. 
----------
Parameters
----------
file1: a .txt file which holds information of our first particles x,y,z 
position for every timestep t. 
file: a .txt file which holds information of our  particles x,y,z 
position for every timestep t. 
file: a .txt file which holds information of our first particles x,y,z 
position for every timestep t. 
file: a .txt file which holds information of our first particles x,y,z 
position for every timestep t. 
file: a .txt file which holds information of our first particles x,y,z 
position for every timestep t. 

"""
def ploterrors(file1,file2,file3,file4,file5, title):
    error, t = errorArray(file1)
    plt.plot(t,error, label=r'h = 0.5')

    error, t = errorArray(file2)
    plt.plot(t,error, label=r'h = 0.1')

    error, t = errorArray(file3)
    plt.plot(t,error, label=r'h = 0.01')

    error, t = errorArray(file4)
    plt.plot(t,error, label=r'h = 0.001')

    error, t = errorArray(file5)
    plt.plot(t,error, label=r'h = 0.0001')

    plt.legend()
    plt.ylim(-0.2,8)
    plt.savefig("../figures/" + title)
    plt.savefig("../figures/" + title)
    plt.show()

"""
Plots the error between our numerical solution and an analytical 
solution with same initial conditions for multiple simulations
with different timesteps. Here we subplot RK4 method vs forward Euleur.
"""
def subplot_errors(file1,file2,file3,file4,file5,file6,file7,file8,file9,file10):

    fig, axs = plt.subplots(nrows=1, ncols=2,figsize=(8,6))


    error, t = errorArray(file1)
    axs[0].plot(t,error, label=r'h = 0.5')

    error, t = errorArray(file2)
    axs[0].plot(t,error, label=r'h = 0.1')

    error, t = errorArray(file3)
    axs[0].plot(t,error, label=r'h = 0.01')

    error, t = errorArray(file4)
    axs[0].plot(t,error, label=r'h = 0.001')

    error, t = errorArray(file5)
    axs[0].plot(t,error, label=r'h = 0.0001')


    error, t = errorArray(file6)
    axs[1].plot(t,error, label=r'h = 0.5')

    error, t = errorArray(file7)
    axs[1].plot(t,error, label=r'h = 0.1')

    error, t = errorArray(file8)
    axs[1].plot(t,error, label=r'h = 0.01')

    error, t = errorArray(file9)
    axs[1].plot(t,error, label=r'h = 0.001')

    error, t = errorArray(file10)
    axs[1].plot(t,error, label=r'h = 0.0001')


    axs[0].set_ylim(0,0.8)
    axs[1].set_ylim(0,0.8)

    fig.text(0.04, 0.5, r'$r_{err}$', ha='center',rotation='vertical')
    #fig.text(0.04, 0.5, r"$v_y$", va='center', rotation='vertical')

    axs[0].set_xlabel(r'time $\mu s$')
    axs[1].set_xlabel(r'time $\mu s$')

    axs[0].legend(loc=1)
    axs[1].legend(loc=1)

    plt.savefig("../figures/RK4_vs_euler.pdf")
    plt.savefig("../figures/RK4_vs_euler.pdf")
    plt.show()

"""
Function to calculate the i-th iteration where our maximum error
occurs. This is to be used for calculating the convergence error.
----------
Parameters
----------
file: a .txt file which holds information of our particles x,y,z 
position for every timestep t. 
"""
def max_error(file):
    x, y, z, t = np.loadtxt(file, unpack=True, usecols=(0,1,2,3))

    B_0 = 1; V_0 = 10; d = 1
    q = 2
    m = 40.07
    x_0 = 1
    z_0 = 1
    v_0 = 0.2
    exX, exY, exZ = analytical3D(t, B_0, V_0, d, q, m, x_0, z_0, v_0)

    test = 0

    #error = np.copy(t)
    max_error = 0
    for i in range(len(t)):
        err_x = np.abs(exX[i] - x[i])
        err_y = np.abs(exY[i] - x[i])
        err_z = np.abs(exZ[i] - x[i])

        #dist_est = np.sqrt(x[i]**2 + y[i]**2 + z[i]**2)
        error = np.abs(np.sqrt(err_x**2 + err_y**2 + err_z**2))
        if error > max_error:
            max_error = error
            test = i
    return test

"""
Function to calculate the convergence error for our numerical method.
----------
Parameters
----------
file1-5: a .txt file which holds information of our particles x,y,z 
position for every timestep t. Each have a unique timestep.
dt: if above 0 we use dt as used for forward euler, if not we use the 
values for RK4. 
"""
def convergence_rate(file1,file2,file3,file4,file5, dt=0):
    max_err_arr = np.zeros(5)
    max_err_arr[0] = max_error(file5)
    max_err_arr[1] = max_error(file4)
    max_err_arr[2] = max_error(file3)
    max_err_arr[3] = max_error(file2)
    max_err_arr[4] = max_error(file1)

    timesteps = np.zeros(5)
    timesteps[0] = 0.1
    timesteps[1] = 0.2
    timesteps[2] = 0.3
    timesteps[3] = 0.4
    timesteps[4] = 0.5

    if (dt > 0):
        timesteps[0] = 0.001
        timesteps[1] = 0.0012
        timesteps[2] = 0.0014
        timesteps[3] = 0.0016
        timesteps[4] = 0.0018


    conv_sum = 0
    for i in range(1,5):
        conv_sum += (np.log(max_err_arr[i]) - np.log(max_err_arr[i-1]))/(np.log(timesteps[i]) - np.log(timesteps[i-1]))

    return conv_sum/4

"""
Function to plot the number of escaped particles for 3 different amplitudes
againt a changing frequency for a total time of 500 mikroseconds.
----------
Parameters
----------
file1-3: a .txt file with the info on number of particles escaped for each 
frequency, for each unique amplitude. 
"""
def plotfreq(file1, file2, file3):
    f1, num_escaped1 = np.loadtxt(file1, unpack=True, usecols=(0,1))
    f2, num_escaped2 = np.loadtxt(file2, unpack=True, usecols=(0,1))
    f3, num_escaped3 = np.loadtxt(file3, unpack=True, usecols=(0,1))
    omega = np.linspace(0.2, 2.5, int(2.3/0.02))

    plt.plot(omega, num_escaped1, label=r'\(f = 0.1\)')
    plt.plot(omega, num_escaped2, label=r'\(f = 0.4\)')
    plt.plot(omega, num_escaped3, label=r'\(f = 0.7\)')
    plt.xlabel(r'\(\omega_V\)')
    plt.xticks(np.arange(0.2, 2.5, 0.25))
    plt.axvline(x=0.63, color='r', linestyle='--', label=r'\(x=0.63\)')
    plt.axvline(x=0.2, color='m', linestyle='--', label=r'\(x=0.2\)')
    plt.axvline(x=0.3, color='c', linestyle='--', label=r'\(x=0.3\)')
    plt.ylabel(r'The number of escaped particles')
    plt.legend()
    plt.show()

"""
Function to plot a zoomed in version of the number of escaped particles 
for 3 different amplitudes againt a changing frequency for a total 
time of 500 mikroseconds.
----------
Parameters
----------
file1-6: a .txt file with the info on number of particles escaped for each 
frequency, for each unique amplitude. 
"""
def plotfreq_zoom(file1, file2, file3, file4, file5, file6):
    f1, num_escaped1 = np.loadtxt(file1, unpack=True, usecols=(0,1))
    f2, num_escaped2 = np.loadtxt(file2, unpack=True, usecols=(0,1))
    f3, num_escaped3 = np.loadtxt(file3, unpack=True, usecols=(0,1))

    f1_, num_escaped1_ = np.loadtxt(file4, unpack=True, usecols=(0,1))
    f2_, num_escaped2_ = np.loadtxt(file5, unpack=True, usecols=(0,1))
    f3_, num_escaped3_ = np.loadtxt(file6, unpack=True, usecols=(0,1))

    omega = np.linspace(0.25, 0.35, int(0.1/0.001))

    fig, ax = plt.subplots(2, 1, sharey=True)

    ax[0].plot(omega, num_escaped1, label=r'\(f = 0.1\)')
    ax[0].plot(omega, num_escaped2, label=r'\(f = 0.4\)')
    ax[0].plot(omega, num_escaped3, label=r'\(f = 0.7\)')

    # ax[1].subplot(1, 2, 2)
    ax[1].plot(omega, num_escaped1_, label=r'\(f = 0.1\)')
    ax[1].plot(omega, num_escaped2_, label=r'\(f = 0.4\)')
    ax[1].plot(omega, num_escaped3_, label=r'\(f = 0.7\)')

    omega = np.linspace(0.25, 0.35, int(0.1/0.001))

    fig.text(0.5, 0.04, r'\(\omega_V\)', ha='center')
    fig.text(0.04, 0.5, r'The number of escaped particles', va='center', rotation='vertical')

    plt.legend()
    plt.show()


"""
Function for solving analytical solution on xy plane
----------
Parameters
----------
t: time we want to solve for.
A_plus: A constant we multiply our exp function with.
A_minus: A constant we multiply our exp function with
w_pluss: Variable dependent on the frequencies w_0 and w_z
w_minus: Variable dependent on the frequencies w_0 and w_z
"""
def analSol2d(t, A_plus, A_minus, w_plus, w_minus):
    position = A_plus*np.exp((-1J)*w_plus*t) + A_minus*np.exp((-1J)*w_minus*t)
    return position.real, position.imag


"""
Function for solving analytical solution on xy plane and plotting it.
----------
Parameters
----------
t: time we want to solve for.
A_plus: A constant we multiply our exp function with.
A_minus: A constant we multiply our exp function with
w_pluss: Variable dependent on the frequencies w_0 and w_z
w_minus: Variable dependent on the frequencies w_0 and w_z
"""
def plotAnalSol2d(A_plus, A_minus, w_plus, w_minus):
    t = np.linspace(0, 100, 10000)
    x = np.zeros(10000)
    y = np.zeros(10000)

    for i in range(len(t)):
        x[i], y[i] = analSol2d(t[i], A_plus, A_minus, w_plus, w_minus)

    plt.plot(x, y, 'c')
    plt.xlabel(r'\(x\)')
    plt.ylabel(r'\(y\)')
    plt.xticks([])
    plt.yticks([])
    plt.show()


"""
Function for solving analytical solution 3 dimentional
----------
Parameters
----------
t: time we want to solve for.
A_plus: A constant we multiply our exp function with.
A_minus: A constant we multiply our exp function with
w_pluss: Variable dependent on the frequencies w_0 and w_z
w_minus: Variable dependent on the frequencies w_0 and w_z
z_0: initial position in z-axis
w_z: frequency for z-axis
"""
def analSol3d(t, A_plus, A_minus, w_plus, w_minus, z_0, w_z):
    position_xy = A_plus*np.exp((-1J)*w_plus*t) + A_minus*np.exp((-1J)*w_minus*t)
    position_z = z_0*np.cos(w_z*t)
    return position_xy.real, position_xy.imag, position_z


"""
Function for solving analytical solution for 3 dimentions and plotting it.
----------
Parameters
----------
t: time we want to solve for.
A_plus: A constant we multiply our exp function with.
A_minus: A constant we multiply our exp function with
w_pluss: Variable dependent on the frequencies w_0 and w_z
w_minus: Variable dependent on the frequencies w_0 and w_z
z_0: initial position in z-axis
w_z: frequency for z-axis
"""
def plotAnalSol3d(A_plus, A_minus, w_plus, w_minus, z_0, w_z):
    t = np.linspace(0, 100, 10000)
    x = np.zeros(10000)
    y = np.zeros(10000)
    z = np.zeros(10000)

    for i in range(len(t)):
        x[i], y[i], z[i] = analSol3d(t[i], A_plus, A_minus, w_plus, w_minus, z_0, w_z)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(x, y, z, c='indigo')
    ax.set_zlim3d([-1.5,2.5])
    ax.grid(False)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    ax.set_xlabel(r'\(x\)')
    ax.set_ylabel(r'\(y\)')
    ax.set_zlabel(r'\(z\)')
    plt.show()
