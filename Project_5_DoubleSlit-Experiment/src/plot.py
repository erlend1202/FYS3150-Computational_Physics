# from typing_extensions import Unpack
import numpy as np
import math
import pyarma as pa
import re

import seaborn as sns

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "sans-serif",
#     "font.sans-serif": ["Helvetica"]})
# # for Palatino and other serif fonts use:
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "serif",
#     "font.serif": ["Palatino"],
# })

"""
Makes a plot of the matrix V. Can be used to check if the slits and borders
are set up correctly
"""
def plot_V():
    h = 0.005
    M = int(1/h) + 1

    filename = "../textfiles/V.bin"
    V = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    V.load(filename)

    V = np.array(V, order="C")

    # plt.imshow(V.real, cmap="magma")
    # plt.show()
    return V



"""
Makes an animation of our simulation.
----------
Parameters
----------
filename: a .bin file containing our values simulated in the experiment. 
T: total time passed in the animation.
"""
def animate(filename, T):
    h = 0.005
    M = int(1/h) + 1

    #
    # Let's generate a dummy time series for a function z(x,y,t)
    #

    # Set up a 2D xy grid
    h = 0.005
    x_points = np.arange(0, 1+h, h)
    y_points = np.arange(0, 1+h, h)
    x, y = np.meshgrid(x_points, y_points, sparse=True)

    # Array of time points
    dt = 2.5*pow(10, -5)
    t_points = np.arange(0, 1+dt, dt)

    # Fill z_data_list with the data_points
    z_data_list = []

    fileloc = "../textfiles/" + filename
    A = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A.load(fileloc)
    A = np.array(A, order="C")
    # print(A)
    timesteps = int(T / (2.5*10**(-5)))
    V = plot_V()
    for i in range(timesteps):
        u = A[:, i].reshape(M-2, M-2)
        P = u.real**2 + u.imag**2
        z_data_list.append(P.T)

    #
    # Now the list z_data_list contains a series of "frames" of z(x,y,t),
    # where each frame can be plotted as a 2D image using imshow. Let's
    # animate it!
    #

    # Some settings
    fontsize = 12
    t_min = t_points[0]
    x_min, x_max = x_points[0], x_points[-1]
    y_min, y_max = y_points[0], y_points[-1]

    # Create figure
    fig = plt.figure()
    ax = plt.gca()
    # print(V.real.shape)
    # Create a colour scale normalization according to the max z value in the first frame
    norm = plt.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[0]))

    # Plot the first frame
    img = ax.imshow(z_data_list[0], extent=[
                    x_min, x_max, y_min, y_max], cmap=plt.get_cmap("magma"), norm=norm)

    # Axis labels
    plt.xlabel("x", fontsize=fontsize)
    plt.ylabel("y", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    # # Add a colourbar
    # cbar = fig.colorbar(img, ax=ax)
    # cbar.set_label("z(x,y,t)", fontsize=fontsize)
    # cbar.ax.tick_params(labelsize=fontsize)

    # Add a text element showing the time
    time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white",
                        horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

    # Function that takes care of updating the z data and other things for each frame

    def animation(i):
        # Normalize the colour scale to the current frame?
        norm = plt.cm.colors.Normalize(
            vmin=0.0, vmax=np.max(np.sqrt(z_data_list[i])))
        img.set_norm(norm)

        # Update z data
        img.set_data(np.sqrt(z_data_list[i]))

        # Update the time label
        current_time = t_min + i * dt
        time_txt.set_text("t = {:.3e}".format(current_time))

        return img

    # Use matplotlib.animation.FuncAnimation to put it all together
    anim = FuncAnimation(fig, animation, interval=20, frames=np.arange(
        0, len(z_data_list), 1), repeat=True, blit=0)

    # Run the animation!
    plt.show()

    # Save the animation
    anim.save('./animation.mp4', writer="ffmpeg", bitrate=-1, fps=30)


"""
Makes an animation of our wave in 3D.
----------
Parameters
----------
filename: a .bin file containing our values simulated in the experiment. 
T: total time passed in the animation.
"""
def animateWave(filename, T):
    h = 0.005
    M = int(1/h) + 1

    #
    # Let's generate a dummy time series for a function z(x,y,t)
    #

    # Set up a 2D xy grid
    h = 0.005
    x = np.linspace(0, 1, M-2)
    x = np.linspace(0, 1, M-2)
    X, Y = np.meshgrid(x, x)

    # Array of time points
    dt = 2.5*pow(10, -5)
    t_points = np.arange(0, 1+dt, dt)

    # Fill z_data_list with the data_points
    z_data_list = []

    fileloc = "../textfiles/" + filename
    A = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A.load(fileloc)
    A = np.array(A, order="C")
    # print(A)
    timesteps = int(T / (2.5*10**(-5)))
    z = np.zeros((M-2, M-2, timesteps))
    for i in range(timesteps):
        u = A[:, i].reshape(M-2, M-2)
        P = u.real**2 + u.imag**2
        z_data_list.append(P.T)
        z[:, :, i] = P.T

    #
    # Now the list z_data_list contains a series of "frames" of z(x,y,t),
    # where each frame can be plotted as a 2D image using imshow. Let's
    # animate it!
    #

    # Create figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    def update_plot(frame_number, zarray, plot):
        plot[0].remove()
        plot[0] = ax.plot_surface(
            x, y, zarray[:, :, frame_number], cmap="magma")

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    N = 14
    nmax = 20
    x = np.linspace(0, 1, M-2)
    x, y = np.meshgrid(x, x)
    zarray = z

    # def f(x, y, sig): return 1/np.sqrt(sig)*np.exp(-(x**2+y**2)/sig**2)

    # for i in range(nmax):
    #     zarray[:, :, i] = f(x, y, 1.5+np.sin(i*2*np.pi/nmax))

    plot = [ax.plot_surface(x, y, zarray[:, :, 0],
                            color='0.75', rstride=1, cstride=1)]
    ax.set_zlim(0, 1.5)
    animate = FuncAnimation(
        fig, update_plot, nmax, fargs=(zarray, plot))
    plt.show()


"""
Makes a subplot containing the normalized probability distribution, and also
a corresponding heatmap for three different simulations.
Preferably the input files should be for a single-slit, double-slit and 
tripple-slit experiment. 
----------
Parameters
----------
filename(1,2,3): a .bin file containing our values simulated in the experiment. 
T: total time passed in the animation.
"""
def detector(filename, filename1, filename2, filename3, T):
    h = 0.005
    M = int(1/h) + 1

    # Set up a 2D xy grid
    h = 0.005
    x_points = np.arange(0, 1+h, h)
    y_points = np.arange(0, 1+h, h)
    x, y = np.meshgrid(x_points, y_points, sparse=True)

    # Array of time points
    dt = 2.5*pow(10, -5)
    t_points = np.arange(0, 1+dt, dt)

    # Fill z_data_list with the data_points
    z_data_list = []
    z_data_list1 = []
    z_data_list2 = []
    z_data_list3 = []

    fileloc = "../textfiles/" + filename
    fileloc1 = "../textfiles/" + filename1
    fileloc2 = "../textfiles/" + filename2
    fileloc3 = "../textfiles/" + filename3

    A = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A.load(fileloc)
    A = np.array(A, order="C")

    A1 = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A1.load(fileloc1)
    A1 = np.array(A1, order="C")

    A2 = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A2.load(fileloc2)
    A2 = np.array(A2, order="C")

    A3 = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A3.load(fileloc3)
    A3 = np.array(A3, order="C")

    timesteps = int(T / (2.5*10**(-5)))
    for i in range(timesteps):
        u = A[:, i].reshape(M-2, M-2)
        P = u.real**2 + u.imag**2
        z_data_list.append(P.T)

        u1 = A1[:, i].reshape(M-2, M-2)
        P1 = u1.real**2 + u1.imag**2
        z_data_list1.append(P1.T)

        u2 = A2[:, i].reshape(M-2, M-2)
        P2 = u2.real**2 + u2.imag**2
        z_data_list2.append(P2.T)

        u3 = A3[:, i].reshape(M-2, M-2)
        P3 = u3.real**2 + u3.imag**2
        z_data_list3.append(P3.T)

    x = np.linspace(0, 1, M-2)

    y_data = z_data_list[-1][:, int(0.8*(M-2))]
    y_data = y_data / sum(y_data)

    y_data1 = z_data_list1[-1][:, int(0.8*(M-2))]
    y_data1 = y_data1 / sum(y_data1)

    y_data2 = z_data_list2[-1][:, int(0.8*(M-2))]
    y_data2 = y_data2 / sum(y_data2)

    y_data3 = z_data_list3[-1][:, int(0.8*(M-2))]
    y_data3 = y_data3 / sum(y_data3)

    # plt.subplot(2, 3, 1)

    mosaic = [['single', 'double'],
              ['single', 'double'],
              ['single_heat', 'double_heat'],
              ['triple', 'triple'],
              ['triple', 'triple'],
              ['triple_heat', 'triple_heat'], ]

    fig = plt.figure(constrained_layout=True)
    ax_dict = fig.subplot_mosaic(mosaic)

    fig.supylabel(r"$p(y\, | \, x = 0.8, \, t=0.002)$")

    # fig.text(0.07, 0.5, r"$p(y\, | \, x = 0.8, \, t=0.002)$",
    #          va='center', rotation='vertical')

    ax_dict['single'].plot(x, y_data1, c="darkmagenta", label="Single-slit")
    ax_dict['single'].legend()

    # ax_dict['single'].plot(x, y_data1, c="darkmagenta", label="Single-Slit")
    # ax_dict['single'].legend()

    # plt.subplot(2, 3, 2)
    ax_dict['double'].plot(x, y_data2, c="midnightblue", label="Double-slit")
    ax_dict['double'].legend()

    # plt.subplot(2, 3, 3)
    ax_dict['triple'].plot(
        x, y_data3, c="mediumslateblue", label="Triple-slit")
    ax_dict['triple'].legend()

    # plt.subplot(2, 3, 4)
    # sns.heatmap([y_data], ax=ax_dict['single_heat'], cbar=False,
    #             xticklabels=False, yticklabels=False)

    sns.heatmap([y_data1], ax=ax_dict['single_heat'], cbar=False,
                xticklabels=False, yticklabels=False)

    # plt.subplot(2, 3, 5)
    sns.heatmap([y_data2], ax=ax_dict['double_heat'], cbar=False,
                xticklabels=False, yticklabels=False)

    # plt.subplot(2, 3, 6)
    sns.heatmap([y_data3], cbar=False, ax=ax_dict['triple_heat'],
                xticklabels=False, yticklabels=False)
    plt.legend([], [], frameon=False)
    # plt.xlabel(r"$x$")
    # axs[1, 0].set(xlabel=r'$x$')
    # axs[1, 1].set(xlabel=r'$x$')
    # axs[1, 2].set(xlabel=r'$x$')
    # plt.savefig("task_9_a_without_log.pdf")

    plt.show()


"""
Reads a .bin file and transforms the array values of u into a matrix. 
----------
Parameters
----------
filename: a .bin file containing our values simulated in the experiment. 
"""
def getMatA(filename):
    # Set up a 2D xy grid
    h = 0.005
    M = int(1/h) + 1

    myfile = "../textfiles/" + filename
    A = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A.load(myfile)
    A = np.array(A, order="C")
    return A


"""
Makes a plot of the complex, real and imaginary part of our simulation, 
for 3 different time points.
----------
Parameters
----------
filename: a .bin file containing our values simulated in the experiment. 
T: total time passed in the animation.
timepoint(2,3): the timepoint we want to plot for.
"""
def plotRealAndImag(filename, T, timePoint, timePoint2, timePoint3):
    h = 0.005
    M = int(1/h) + 1
    dt = 2.5*pow(10, -5)

    A = getMatA(filename)

    # Fill z_data_list with the data_points
    z_data_list = []

    # Stores lists for real and imaginary data points
    real_data_list = []
    imag_data_list = []

    # print(A)
    # Getting all frames from simulation
    timesteps = int(T / (2.5*10**(-5)))
    for i in range(timesteps):
        u = A[:, i].reshape(M-2, M-2)
        P = u.real**2 + u.imag**2
        R = u.real
        I = u.imag

        z_data_list.append(P.T)
        real_data_list.append(R.T)
        imag_data_list.append(I.T)

    fig, axs = plt.subplots(3, 3, figsize=(8, 8))

    plt.rcParams.update({'font.size': 8})
    # Task 8
    # ----
    if timePoint == 0.002:
        timestamp = 78
    else:
        timestamp = int(timePoint/dt)
    # -----

    timestamp = [int(timePoint//dt), int(timePoint2//dt), int(timePoint3//dt)]

    sns.heatmap(z_data_list[timestamp[0]], cbar=False, ax=axs[0, 0],
                xticklabels=False, yticklabels=False)

    sns.heatmap(real_data_list[timestamp[0]], cbar=False, ax=axs[1, 0],
                xticklabels=False, yticklabels=False)

    sns.heatmap(imag_data_list[timestamp[0]], cbar=False, ax=axs[2, 0],
                xticklabels=False, yticklabels=False)

    sns.heatmap(z_data_list[timestamp[1]], cbar=False, ax=axs[0, 1],
                xticklabels=False, yticklabels=False)

    sns.heatmap(real_data_list[timestamp[1]], cbar=False, ax=axs[1, 1],
                xticklabels=False, yticklabels=False)

    sns.heatmap(imag_data_list[timestamp[1]], cbar=False, ax=axs[2, 1],
                xticklabels=False, yticklabels=False)

    sns.heatmap(z_data_list[timestamp[2]], cbar=False, ax=axs[0, 2],
                xticklabels=False, yticklabels=False)

    sns.heatmap(real_data_list[timestamp[2]], cbar=False, ax=axs[1, 2],
                xticklabels=False, yticklabels=False)

    sns.heatmap(imag_data_list[timestamp[2]], cbar=False, ax=axs[2, 2],
                xticklabels=False, yticklabels=False)

    # fig.text(0.5, 0.2, 'x', ha='center', va='center', size='large')
    # fig.text(0.01, 0.5, r'$y$', ha='center', va='center',
    #          rotation='vertical', size='large')
    plt.tight_layout()
    fig.savefig("../figures/task8_.png")
    plt.show()

"""
Makes a plot of the development probability, with and without a double-slit
barrier.
----------
Parameters
----------
filename: a .bin file containing our values simulated in the experiment. 
T: total time passed in the animation.
timepoint(2): the timepoint we want to plot for.
"""
def plotProbDev(filename, filename2, T):
    h = 0.005
    M = int(1/h) + 1

    # Set up a 2D xy grid
    h = 0.005
    x_points = np.arange(0, 1+h, h)
    y_points = np.arange(0, 1+h, h)
    x, y = np.meshgrid(x_points, y_points, sparse=True)

    # Array of time points
    dt = 2.5*pow(10, -5)
    t_points = np.arange(0, 1+dt, dt)

    # Fill z_data_list with the data_points
    z_data_list = []

    fileloc = "../textfiles/" + filename
    fileloc2 = "../textfiles/" + filename2

    A = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A.load(fileloc)
    A = np.array(A, order="C")

    A2 = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A2.load(fileloc2)
    A2 = np.array(A2, order="C")

    timesteps = int(T / (2.5*10**(-5)))
    # Array for storting total prob at all timesteps
    totProb = np.zeros(timesteps)
    totProb2 = np.zeros(timesteps)

    for i in range(timesteps):
        u = A[:, i].reshape(M-2, M-2)
        P = u.real**2 + u.imag**2
        totProb[i] = 1 - np.sum(P)

        u2 = A2[:, i].reshape(M-2, M-2)
        P2 = u2.real**2 + u2.imag**2
        totProb2[i] = 1 - np.sum(P2)

    x = np.linspace(0, T, timesteps)

    plt.rcParams.update({'font.size': 14})

    plt.subplot(1, 2, 1)
    plt.scatter(x, abs(totProb), color="darkmagenta",
                marker=".", label="No Barrier")
    plt.ylim([-0.5*10**(-15), 6*10**(-15)])
    plt.ylabel("Deviation of the Total Probability from 1.0")
    plt.xlabel("Time")
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.scatter(x, abs(totProb2), color="midnightblue",
                marker=".", label="Double-Slit")
    plt.ylim([-0.5*10**(-15), 6*10**(-15)])
    plt.xlabel("Time")
    plt.legend()

    plt.show()

"""
testing 
"""
def plotWaves(filename, T, timePoint):
    h = 0.005
    M = int(1/h) + 1

    # Set up a 2D xy grid
    h = 0.005

    # Array of time points
    dt = 2.5*pow(10, -5)
    t_points = np.arange(0, 1+dt, dt)

    # Fill z_data_list with the data_points
    z_data_list = []

    # Stores lists for real and imaginary data points
    real_data_list = []
    imag_data_list = []

    myfile = "../textfiles/" + filename
    A = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A.load(myfile)
    A = np.array(A, order="C")

    # print(A)
    # Getting all frames from simulation
    timesteps = int(T / (2.5*10**(-5)))
    for i in range(timesteps):
        u = A[:, i].reshape(M-2, M-2)
        P = u.real**2 + u.imag**2
        R = u.real
        I = u.imag

        z_data_list.append(np.sqrt(P.T))
        real_data_list.append(R.T)
        imag_data_list.append(I.T)

    plt.rcParams.update({'font.size': 8})

    x_points = np.linspace(0, 1, M-2)
    y_points = np.linspace(0, 1, M-2)
    X, Y = np.meshgrid(x_points, y_points)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X, Y, P, cmap="magma")
    plt.show()

"""
testing
"""
def hitWall(filename, filename2, filename3, T):
    h = 0.005
    M = int(1/h) + 1

    # Fill z_data_list with the data_points
    z_data_list = []
    z_data_list2 = []
    z_data_list3 = []

    fileloc = "../textfiles/" + filename
    fileloc2 = "../textfiles/" + filename2
    fileloc3 = "../textfiles/" + filename3

    A = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A.load(fileloc)
    A = np.array(A, order="C")

    A2 = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A2.load(fileloc2)
    A2 = np.array(A2, order="C")

    A3 = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A3.load(fileloc3)
    A3 = np.array(A3, order="C")

    timesteps = int(T / (2.5*10**(-5)))
    for i in range(timesteps):
        u = A[:, i].reshape(M-2, M-2)
        P = u.real**2 + u.imag**2
        z_data_list.append(P.T)

        u2 = A2[:, i].reshape(M-2, M-2)
        P2 = u2.real**2 + u2.imag**2
        z_data_list2.append(P2.T)

        u3 = A3[:, i].reshape(M-2, M-2)
        P3 = u3.real**2 + u3.imag**2
        z_data_list3.append(P3.T)

    x = np.linspace(0, 1, M-2)

    y_data = z_data_list[-1][:, int((M-3))]
    y_data = y_data / sum(y_data)

    y_data2 = z_data_list2[-1][:, int((M-3))]
    y_data2 = y_data2 / sum(y_data2)

    y_data3 = z_data_list3[-1][:, int((M-3))]
    y_data3 = y_data3 / sum(y_data3)

    # print(z_data_list[-1][:, int(0.8*(M-2))])

    fig, axs = plt.subplots(3, 1, figsize=(12, 6))
    fig.text(0.05, 0.6, r"$p(y\, | \, x = 0.8, \, t=0.002)$",
             va='center', rotation='vertical')

    # plt.subplot(2, 3, 4)
    sns.heatmap([y_data], ax=axs[0], cbar=False,
                xticklabels=False, yticklabels=False)

    # plt.subplot(2, 3, 5)
    sns.heatmap([y_data2], ax=axs[1], cbar=False,
                xticklabels=False, yticklabels=False)

    # plt.subplot(2, 3, 6)
    sns.heatmap([y_data3], cbar=False, ax=axs[2],
                xticklabels=False, yticklabels=False)
    plt.legend([], [], frameon=False)
    # plt.xlabel(r"$x$")
    # axs[1, 0].set(xlabel=r'$x$')
    # axs[1, 1].set(xlabel=r'$x$')
    # axs[1, 2].set(xlabel=r'$x$')
    # plt.savefig("task_9_a_without_log.pdf")

    plt.show()

"""
testing
"""
def plotMultSlits(filename, filename1, filename2, filename3, T, timePoint, timePoint2, timePoint3):
    h = 0.005
    M = int(1/h) + 1

    # Set up a 2D xy grid
    h = 0.005
    x_points = np.arange(0, 1+h, h)
    y_points = np.arange(0, 1+h, h)
    x, y = np.meshgrid(x_points, y_points, sparse=True)

    # Array of time points
    dt = 2.5*pow(10, -5)
    t_points = np.arange(0, 1+dt, dt)

    # Fill z_data_list with the data_points
    z_data_list = []
    z_data_list1 = []
    z_data_list2 = []
    z_data_list3 = []

    fileloc = "../textfiles/" + filename
    fileloc1 = "../textfiles/" + filename1
    fileloc2 = "../textfiles/" + filename2
    fileloc3 = "../textfiles/" + filename3

    A = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A.load(fileloc)
    A = np.array(A, order="C")

    A1 = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A1.load(fileloc1)
    A1 = np.array(A1, order="C")

    A2 = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A2.load(fileloc2)
    A2 = np.array(A2, order="C")

    A3 = pa.cx_mat()  # Create pa.mat object (just as arma::mat in C++)
    # Load the content of the matrix you saved into your Python program.
    A3.load(fileloc3)
    A3 = np.array(A3, order="C")

    timesteps = int(T / (2.5*10**(-5)))
    for i in range(timesteps):
        u = A[:, i].reshape(M-2, M-2)
        P = u.real**2 + u.imag**2
        z_data_list.append(P.T)

        u1 = A1[:, i].reshape(M-2, M-2)
        P1 = u1.real**2 + u1.imag**2
        z_data_list1.append(P1.T)

        u2 = A2[:, i].reshape(M-2, M-2)
        P2 = u2.real**2 + u2.imag**2
        z_data_list2.append(P2.T)

        u3 = A3[:, i].reshape(M-2, M-2)
        P3 = u3.real**2 + u3.imag**2
        z_data_list3.append(P3.T)

    fig, axs = plt.subplots(3, 3, figsize=(8, 8))

    timestamp = [int(timePoint//dt), int(timePoint2//dt), int(timePoint3//dt)]

    sns.heatmap(z_data_list1[timestamp[0]], cbar=False, ax=axs[0, 0],
                xticklabels=False, yticklabels=False)

    sns.heatmap(z_data_list2[timestamp[0]], cbar=False, ax=axs[1, 0],
                xticklabels=False, yticklabels=False)

    sns.heatmap(z_data_list3[timestamp[0]], cbar=False, ax=axs[2, 0],
                xticklabels=False, yticklabels=False)

    sns.heatmap(z_data_list1[timestamp[1]], cbar=False, ax=axs[0, 1],
                xticklabels=False, yticklabels=False)

    sns.heatmap(z_data_list2[timestamp[1]], cbar=False, ax=axs[1, 1],
                xticklabels=False, yticklabels=False)

    sns.heatmap(z_data_list3[timestamp[1]], cbar=False, ax=axs[2, 1],
                xticklabels=False, yticklabels=False)

    sns.heatmap(z_data_list1[timestamp[2]], cbar=False, ax=axs[0, 2],
                xticklabels=False, yticklabels=False)

    sns.heatmap(z_data_list2[timestamp[2]], cbar=False, ax=axs[1, 2],
                xticklabels=False, yticklabels=False)

    sns.heatmap(z_data_list3[timestamp[2]], cbar=False, ax=axs[2, 2],
                xticklabels=False, yticklabels=False)

    # fig.text(0.5, 0.2, 'x', ha='center', va='center', size='large')
    # fig.text(0.01, 0.5, r'$y$', ha='center', va='center',
    #          rotation='vertical', size='large')
    plt.tight_layout()
    # fig.savefig("../figures/task8_.png")
    plt.show()

""" 
Calls plotProbDev.
----------
Parameters
----------
filename: a .bin file containing our values simulated in the experiment. 
"""
def oppgave7(filenames):
    plotProbDev(*filenames, 0.008)

""" 
Calls plotRealAndImag.
----------
Parameters
----------
filename: a .bin file containing our values simulated in the experiment. 
timepoints: list of timepoints we want to use.
"""
def oppgave8(filename, timepoints):
    plotRealAndImag(filename, 0.002, *timepoints)

""" 
Calls detector().
----------
Parameters
----------
filename: a .bin file containing our values simulated in the experiment. 
"""
def oppgave9(filenames):
    detector(*filenames, 0.002)

""" 
testing 
"""
def waves(filename, timepoints):
    plotWaves(filename, 0.002, timepoints[0])


if __name__ == '__main__':
    # noSlit = "results_no_slit.bin"
    singleSlit = "results_single_slit.bin"
    doubleSlit = "results_double_slit.bin"
    # doubleSlit_task8 = "results_double_slit_task_8.bin"
    tripleSlit = "results_triple_slit.bin"
    # doubleSlit_animation = "results_double_slit_animation.bin"

    # problem7Configs = [noSlit, doubleSlit]
    # oppgave7(problem7Configs)

    # Husk å lage ny double slit me sigma_y = 0.2
    # plt.rcParams.update({'font.size': 8})
    problem8Configs = [0, 0.001, 0.002]
    # oppgave8(doubleSlit, problem8Configs)
    # waves(doubleSlit, problem8Configs)
    # plotMultConfigs = [noSlit, singleSlit, doubleSlit, tripleSlit]
    # plotMultSlits(*plotMultConfigs, 0.002, *problem8Configs)

    # plt.rcParams.update({'font.size': 13})
    problem9Configs = [singleSlit, singleSlit, doubleSlit, tripleSlit]
    oppgave9(problem9Configs)
    plot_V()

    animate(doubleSlit_animation, 0.015)
    # detector(doubleSlit_task8, 0.002)

    # hitWall(*problem9Configs, 0.0028)

    # # # #For task 8
    # plotRealAndImag(singleSlit, 0.008, 0.0)
    # plotRealAndImag(singleSlit, 0.008, 0.001)
    # plotRealAndImag(singleSlit, 0.008, 0.002)
    # plotProbDev(noSlit, 0.008)
