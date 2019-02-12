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

def draw_fracture(simulation_id):

    plt.figure(simulation_id)
    plt.scatter(50*np.sqrt(3), 0.005, marker="+")

#------------------------------------------------------------------------------#
# add data to plots

# add plots for the different schemes
# The first argument is the path to your data, where the data is assumed to be ordered comma-separated in the following order:
#   -> 1. arc length, 2. value of either c or pressure
# The second argument defines the legend you want to add to your data (Institution / Numerical method)

titles = ['$\\sim 30k$ cells', '$\\sim 150k$ cells']
refinement_index = ['0', '1']

for title, ref in zip(titles, refinement_index):

    # University of Stuttgart Mpfa results
    # TODO add coarse solution data
    if ref == str(0):
        fileNameBody = "../results/USTUTT/MPFA/dol_line_"
        plot(fileNameBody+"0_refinement_"+ref+".csv", "USTUTT-MPFA", id_p_0_matrix, title)
        plot(fileNameBody+"1_refinement_"+ref+".csv", "USTUTT-MPFA", id_p_1_matrix, title)

    # University of Bergen
    methods = ['MVEM', 'TPFA', 'MPFA', 'RT0']
    for method in methods:
        folder = "../results/UiB/"+method+"/"
        label = "UiB-"+method
        plot(folder+"dol_line_0_refinement_"+ref+".csv", label, id_p_0_matrix, title)
        plot(folder+"dol_line_1_refinement_"+ref+".csv", label, id_p_1_matrix, title)

    # UNICE-UNIGE
    methods = ['HFVCont', 'HFVDisc', 'VAGCont', 'VAGDisc']
    for method in methods:
        folder = "../results/UNICE-UNIGE/"+method+"/"
        label = "UNICE-UNIGE-"+method
        plot(folder+"dol_line_0_refinement_"+ref+".csv", label, id_p_0_matrix, title)
        plot(folder+"dol_line_1_refinement_"+ref+".csv", label, id_p_1_matrix, title)

    #draw_fracture(id_c_matrix)

    # save figures
    save(id_p_0_matrix, "pol_p_0_matrix_"+ref)
    save(id_p_1_matrix, "pol_p_1_matrix_"+ref)

#------------------------------------------------------------------------------#
