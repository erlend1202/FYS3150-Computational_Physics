import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt


from scipy.interpolate import UnivariateSpline
from scipy import stats

def plot_burn_in_time(save=False):
    """Reads the output from the function compare_burn_in_time() from main.cpp and
    plots the figures used in the report. The ordered state is anti ferromagnetic (maby try with ferromagnetic?)

    Args:
        save (bool, optional): Saves the figures if set to True. Default set to False.
    """

    rand_file1 = "../textfiles/T1_mean_energy_and_abs_m_random.txt"
    rand_file2 = "../textfiles/T2_mean_energy_and_abs_m_random.txt"

    ord_file1 = "../textfiles/T1_mean_energy_and_abs_m_ordered.txt"
    ord_file2 = "../textfiles/T2_mean_energy_and_abs_m_ordered.txt"

    rand_eps1, rand_m1, cycles = np.loadtxt(rand_file1, unpack=True, usecols=(0,1,2))
    rand_eps2, rand_m2 = np.loadtxt(rand_file2, unpack=True, usecols=(0,1))

    ord_eps1, ord_m1 = np.loadtxt(ord_file1, unpack=True, usecols=(0,1))
    ord_eps2, ord_m2 = np.loadtxt(ord_file2, unpack=True, usecols=(0,1))

    N = 20*20

    rand_eps1, rand_m1, rand_eps2, rand_m2 = rand_eps1, rand_m1, rand_eps2, rand_m2
    ord_eps1, ord_m1, ord_eps2, ord_m2 = ord_eps1, ord_m1, ord_eps2, ord_m2

    plt.figure(figsize=(8,6))
    plt.rcParams.update({'font.size': 14})

    plt.subplot(211)
    plt.plot(cycles, rand_eps1/N, c="C0", label=r"$T = 1 \; J / k_B$, unordered.")
    # plt.scatter(cycles, rand_eps1, c="C0", marker=".")
    plt.plot(cycles, ord_eps1/N, c="C4", label=r"$T = 1 \; J / k_B$, ordered.")
    # plt.scatter(cycles, ord_eps1, c="C4", marker=".")

    plt.plot(cycles, rand_eps2/N, c="C6", label=r"$T = 2.4 \; J / k_B$, unordered.")
    # plt.scatter(cycles, rand_eps2, c="C6", marker=".")
    plt.plot(cycles, ord_eps2/N, c="C9", label=r"$T = 2.4 \; J / k_B$, ordered.")
    # plt.scatter(cycles, ord_eps2, c="C9", marker=".")
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    plt.ylabel(r"$\langle\epsilon\rangle \, \,$ [J]")
    plt.xscale("log")
    plt.legend(fontsize=12, loc="upper right")

    plt.subplot(212)
    plt.plot(cycles, rand_m1/N, c="C0", label=r"$T = 1 \; J / k_B$, unordered.")
    # plt.scatter(cycles, rand_m1, c="C0", marker=".")
    plt.plot(cycles, ord_m1/N, c="C4", label=r"$T = 1 \; J / k_B$, ordered.")
    # plt.scatter(cycles, ord_m1, c="C4", marker=".")


    plt.plot(cycles, rand_m2/N, c="C6", label=r"$T = 2.4 \; J / k_B$, unordered.")
    # plt.scatter(cycles, rand_m2, c="C6", marker=".")
    plt.plot(cycles, ord_m2/N, c="C9", label=r"$T = 2.4 \; J / k_B$, ordered.")
    # plt.scatter(cycles, ord_m2, c="C9", marker=".")
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    plt.xlabel("The number of MC-cycles")
    plt.ylabel(r"$\langle|m|\rangle$")
    plt.legend(fontsize=12)
    plt.xscale("log")

    if save:
        plt.savefig("../figures/burn_in_times.png")

    plt.show()


def plot_hist_pdf(save=False):
    """Read a file containing approximated values of p_eps(eps;T) for T = 1 and T = 2.4. The lattice size is 20.

    Args:
        save (bool, optional): Saves the figures if set to True. Default set to False.
    """
    eps1, eps2 = np.loadtxt("../textfiles/problem6_test.txt", unpack=True, usecols=(0,1))

    fig, axs = plt.subplots(2)
    axs[0].hist(eps1, density=True, bins = 11, ec="white",color='purple', label=r"$T = 1 \; J / k_B$")
    axs[1].hist(eps2, density=True, bins = 54, ec="white",color='orange', label=r"$T = 2.4 \; J/ k_B$")


    axs[0].legend(fontsize=14, loc = 1)
    axs[1].legend(fontsize=14, loc = 1)


    plt.rcParams.update({'font.size': 14})
    fig.text(0.5, 0.04, r"$\epsilon$ [J]", ha='center')
    fig.text(0.04, 0.5, 'Probability [\%]', va='center', rotation='vertical')

    if save:
        plt.savefig("../figures/energy_per_spin_hist.pdf")

    var1 = np.var(eps1)
    var2 = np.var(eps2)

    print(var1)
    print(var2)

    print(var2/var1)

    N = 20*20
    t = np.linspace(0, 10**6, 10**6)
    anal = analytical_energy(t) / N


    plt.show()



def plot_lat_comp(save=False):
    """Reads the files created by laticeSizeComp() and plots
    the data along the temperature range
    """

    # Arrays for all values
    eps = [[] for _ in range(4)]
    mag = [[] for _ in range(4)]
    heat = [[] for _ in range(4)]
    susc = [[] for _ in range(4)]
    temp = [[] for _ in range(4)]

    latices = [40, 60, 80, 100] # Lattice sizes of the given data

    # Read data into arrays
    for i in range(len(latices)):
        filename = f"../textfiles/mat{latices[i]}info.txt"
        eps[i], mag[i], heat[i], susc[i], temp[i] = np.loadtxt(filename, unpack=True, usecols=(0,1,2,3,4))

    for i in range(len(latices)):
        heat[i] = heat[i]


    critical_temp = []  # Array for storing critical temp for different lattice sizes


    # Plotting datapoints
    plt.rcParams.update({'font.size': 14})

    plt.figure(figsize=(8,6))
    plt.subplot(121)
    for i in range(len(latices)):
        plt.scatter(temp[i], eps[i], label=fr"$L={latices[i]}$")

    plt.ylabel(r"$\langle\epsilon\rangle \, \, [J]$")
    plt.xlabel(r"$T \, \, [J / k_B]$")
    plt.legend()

    plt.subplot(122)
    for i in range(len(latices)):
        N = latices[i]**2
        plt.scatter(temp[i], mag[i], label=fr"$L={latices[i]}$")

    plt.ylabel(r"$\langle|m|\rangle$")
    plt.xlabel(r"$T \, \, [J / k_B]$")
    plt.legend()
    plt.tight_layout()

    plt.figure(figsize=(8,6))
    plt.subplot(121)

    smoothings = [.01, .01, 0.01, 0.01] # Smoothing parameters for splines


    # Calculating splines for heat data on for all lattice sizes
    for i in range(len(latices)):
        plt.scatter(temp[i], heat[i], label=fr"$L={latices[i]}$", alpha=0.4)

        x_line = np.linspace((temp[0][0]), (temp[0][-1]), 150000)
        spline = UnivariateSpline(temp[i], heat[i], s=smoothings[i])

        plt.plot(x_line, spline(x_line))

        y = spline(x_line)

        plt.scatter(x_line[np.argmax(y)], y[np.argmax(y)], marker="x", color="k")

        critical_temp.append(x_line[int(np.argmax(y))])



    # Setting labels
    plt.ylabel(r"$C_V \, \, [k_B]$")
    plt.xlabel(r"$T \, \, [J / k_B]$")
    plt.legend()

    smoothings = [10, 50, 200, 250] # Smoothing parameters for splines
    
    plt.subplot(122)

    # Calculating and plotting splines for all lattice sizes
    for i in range(len(latices)):
        plt.scatter(temp[i], susc[i], label=fr"$L={latices[i]}$", alpha=0.4)
        x_line = np.linspace((temp[0][0]), (temp[0][-1]), 150000)
        spline = UnivariateSpline(temp[i], susc[i], s=smoothings[i], k=2)
        plt.plot(x_line, spline(x_line))

        y = spline(x_line)
        plt.scatter(x_line[np.argmax(y)], y[np.argmax(y)], marker="x", color="k")

        critical_temp[i] += x_line[int(np.argmax(y))]
        critical_temp[i] /= 2



    # Setting labels
    plt.ylabel(r"$\chi \, \, [1/J]$")
    plt.xlabel(r"$T \, \, [J / k_B]$")
    plt.legend()
    plt.tight_layout()

    plt.show()


    # Plotting critical temperatures
    inv_latices = np.array(latices)
    inv_latices = 1/inv_latices

    critical_temp = np.array(critical_temp)
    for i, (inv, t) in enumerate(zip(inv_latices, critical_temp)):
        plt.scatter(inv, t, label=fr"$L = ${latices[i]}")


    # Calculating critical temp for inifinte lattice
    res = stats.linregress(inv_latices, critical_temp)
    x_line = np.linspace(-0.005,inv_latices[0]+ 0.005)


    # Setting labels and plotting
    plt.plot(x_line, res.intercept + res.slope*x_line, 'r', label='Linear regression')
    val_test = res.intercept + res.slope*x_line
    plt.plot(0, res.intercept, 'k+', markersize=18, label=fr"$L = \infty$")
    plt.xlabel(r"$L^{-1}$")
    plt.ylabel(r"$T \, \, [J/k_B]$")
    plt.legend()
    plt.show()

    print("Actual critical temperature = ",res.intercept)




def analytical_energy(T):
    beta = 1/T
    val = -8*np.sinh(beta*8)/(np.cosh(beta*8) + 3)
    return val

def analytical_energy2(T):
    beta = 1/T
    val = 64*np.cosh(beta*8)/(np.cosh(beta*8) + 3)
    return val

def analytical_magnetism(T):
    beta = 1/T
    val = (2*np.exp(beta*8) + 4)/(np.cosh(8*beta) + 3)
    return val

def analytical_magnetism2(T):
    beta = 1/T
    val = 8*(np.exp(beta*8) + 1)/(np.cosh(8*beta) + 3)
    return val

def analytical_C(T, E2, E):
    val = 1/(T**2)*( E2 - E*E)
    return val

def analytical_X(T, M2, abs_M):
    val = (1/T)*(M2 - abs_M*abs_M)
    return val


def plot_error(save=False):
    E, E2, T = np.loadtxt("../textfiles/numerical_energy_T.txt", unpack=True, usecols=(0,1,2))
    abs_M, M2 = np.loadtxt("../textfiles/numerical_magnetism_T.txt", unpack=True, usecols=(0,1))
    specific_heat, susceptibility = np.loadtxt("../textfiles/numerical_spec_susc_T.txt", unpack=True, usecols=(0,1))

    plt.figure(figsize=(8,6))
    plt.rcParams.update({'font.size': 14})

    N=4

    plt.subplot(221)
    plt.plot(T, analytical_energy(T), c="darkgreen", label="Analytical")
    plt.scatter(T, E, c="crimson", alpha=0.7, label="Numerical")
    plt.ylabel(r"$\langle E\rangle$ [ J ]")


    plt.subplot(222)
    plt.plot(T, analytical_energy2(T), c="darkgreen", label="Analytical")
    plt.scatter(T, E2, c="crimson", alpha=0.7, label="Numerical")
    plt.ylabel(r"$\langle E^2\rangle$ [ $J^2$ ]")

    plt.subplot(223)
    plt.plot(T, analytical_magnetism(T), c="darkgreen", label="Analytical")
    plt.scatter(T, abs_M, c="crimson", alpha=0.7, label="Numerical")
    plt.ylabel(r"$\langle |M|\rangle$")
    plt.xlabel(r"T [ J / $k_B$ ]")

    plt.subplot(224)
    plt.plot(T, analytical_magnetism2(T), c="darkgreen", label="Analytical")
    plt.scatter(T, M2, c="crimson", alpha=0.7, label="Numerical")
    plt.ylabel(r"$\langle M^2\rangle$")
    plt.xlabel(r"T [ J / $k_B$ ]")

    plt.legend()
    plt.tight_layout()


    if save:
        plt.savefig("../figures/error_compare.pdf")

    plt.show()

def plot_more_error(save=False):
    E, M, cycles = np.loadtxt("../textfiles/error_analysis.txt", unpack=True, usecols=(0,1,2))

    plt.figure(figsize=(8,6))
    plt.rcParams.update({'font.size': 22})

    plt.plot(cycles, np.abs(E-analytical_energy(1))/np.abs(E), "b", label=r"$\langle E\rangle$")
    plt.plot(cycles, np.abs(M-analytical_magnetism(1))/np.abs(M), "r", label=r"$\langle M \rangle$")
    plt.xscale("log")
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel(r"Number of MC cycles")
    plt.ylabel(r"Relative error")
    plt.legend()
    plt.tight_layout()

    if save:
        plt.savefig("../figures/relative_error.pdf")

    plt.show()






if __name__ == '__main__':
    pass