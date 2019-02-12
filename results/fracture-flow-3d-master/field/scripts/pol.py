from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np
import os
import sys

#------------------------------------------------------------------------------#

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=15)

# ids of the different plots
id_p_0_matrix = 0   # pressure along (0, 100, 100)-(100, 0, 0)
id_p_1_matrix = 1   # p along (0, 100, 100)-(100, 0, 0)

def plot(file_name, legend, ID, title):

    c = lambda s: float(s.decode().replace('D', 'e'))
    N = 2
    data = np.genfromtxt(file_name, delimiter=",", converters=dict(zip(range(N), [c]*N)))

    plt.figure(ID)
    plt.plot(data[:, 0], data[:, 1], label=legend)
    plt.xlabel('arc length')
    plt.grid(True)
    plt.title(title)
    plt.legend()
    plt.ylabel(r"$p_3$")

def save(simulation_id, filename):
    folder = "./plots/"

    if not os.path.exists(folder):
        os.makedirs(folder)

    plt.figure(simulation_id)
    plt.savefig(folder+filename, bbox_inches='tight')
    plt.gcf().clear()

#------------------------------------------------------------------------------#

# add data to plots

# add plots for the different schemes
# The first argument is the path to your data, where the data is assumed to be ordered comma-separated in the following order:
#   -> 1. arc length, 2. value of either c or pressure
# The second argument defines the legend you want to add to your data (Institution / Numerical method)

# University of Bergen
methods = ['MVEM', 'TPFA', 'MPFA', 'RT0']
for method in methods:
    folder = "../results/UiB/"+method+"/"
    label = "UiB-"+method
    title = "pressure over line 0"
    plot(folder+"dol_line_0.csv", label, id_p_0_matrix, title)
    title = "pressure over line 1"
    plot(folder+"dol_line_1.csv", label, id_p_1_matrix, title)

# UNICE-UNIGE
methods = ['HFVCont', 'HFVDisc', 'VAGCont', 'VAGDisc']
for method in methods:
    folder = "../results/UNICE-UNIGE/"+method+"/"
    label = "UNICE-UNIGE-"+method
    title = "pressure over line 0"
    plot(folder+"dol_line_0.csv", label, id_p_0_matrix, title)
    title = "pressure over line 1"
    plot(folder+"dol_line_1.csv", label, id_p_1_matrix, title)

    #draw_fracture(id_c_matrix)

# save figures
save(id_p_0_matrix, "pol_line_0")
save(id_p_1_matrix, "pol_line_1")

#------------------------------------------------------------------------------#
