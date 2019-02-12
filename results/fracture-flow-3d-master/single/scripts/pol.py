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
id_p_matrix = 0   # pressure along (0, 100, 100)-(100, 0, 0)
id_c_matrix = 1   # c along (0, 100, 100)-(100, 0, 0)
id_c_fracture = 2 # c along (0, 100, 80)-(100, 0, 20)

def plot(file_name, legend, ID, title):

    c = lambda s: float(s.decode().replace('D', 'e'))
    N = 5
    data = np.genfromtxt(file_name, delimiter=",", converters=dict(zip(range(N), [c]*N)))

    plt.figure(ID)
    plt.plot(data[:, 2*ID], data[:, 2*ID+1], label=legend)
    plt.xlabel('arc length')
    plt.grid(True)
    plt.title(title)
    if ID == id_p_matrix:
        plt.legend(loc="upper right")
    elif ID == id_c_fracture:
        plt.legend(loc="lower left")
    else:
        plt.legend(loc="lower right")

    # choose y-label depending on plot id
    if ID == id_p_matrix:
        plt.ylabel(r"$p_3$")
    elif ID == id_c_matrix:
        plt.ylabel(r"$c_3$")
    elif ID == id_c_fracture:
        plt.ylabel(r"$c_2$")
    else:
        print("Error. Invalid plot id provided.")
        sys.exit(1)

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

titles = ['$\\sim 1k$ cells', '$\\sim 10k$ cells', '$\\sim 100k$ cells']
refinement_index = ['0', '1', '2']

for title, ref in zip(titles, refinement_index):

    # TODO: add reference solution to plots as soon as available
    # University of Stuttgart Mpfa results
    folder = "../results/USTUTT/MPFA/"
    plot(folder+"dol_refinement_"+ref+".csv", "USTUTT-MPFA", id_p_matrix, title)
    plot(folder+"dol_refinement_"+ref+".csv", "USTUTT-MPFA", id_c_matrix, title)
    plot(folder+"dol_refinement_"+ref+".csv", "USTUTT-MPFA", id_c_fracture, title)

    # University of Bergen
    methods = ['MVEM', 'TPFA', 'MPFA', 'RT0']
    for method in methods:
        folder = "../results/UiB/"+method+"/"
        label = "UiB-"+method
        plot(folder+"dol_refinement_"+ref+".csv", label, id_p_matrix, title)
        plot(folder+"dol_refinement_"+ref+".csv", label, id_c_matrix, title)
        plot(folder+"dol_refinement_"+ref+".csv", label, id_c_fracture, title)

    # UNICE-UNIGE
    methods = ['HFVCont', 'HFVDisc', 'VAGCont', 'VAGDisc']
    for method in methods:
        folder = "../results/UNICE-UNIGE/"+method+"/"
        label = "UNICE-UNIGE-"+method
        plot(folder+"dol_refinement_"+ref+".csv", label, id_p_matrix, title)
        plot(folder+"dol_refinement_"+ref+".csv", label, id_c_matrix, title)
        plot(folder+"dol_refinement_"+ref+".csv", label, id_c_fracture, title)

    # NCU_TW
    methods = ["NAME"]
    for method in methods:
        folder = "../results/NCU_TW/"#+method+"/"
        label = "NCU-TV-"+method
        plot(folder+"dol_refinement_"+ref+".csv", label, id_p_matrix, title)
        plot(folder+"dol_refinement_"+ref+".csv", label, id_c_matrix, title)
        plot(folder+"dol_refinement_"+ref+".csv", label, id_c_fracture, title)

    # save figures
    save(id_p_matrix, "pol_p_matrix_"+ref)
    save(id_c_matrix, "pol_c_matrix_"+ref)
    save(id_c_fracture, "pol_c_fracture_"+ref)

#------------------------------------------------------------------------------#
