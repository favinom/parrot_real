import matplotlib.pyplot as plt
import numpy as np
import os
import sys

#------------------------------------------------------------------------------#

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=15)

def plot(file_name, legend, title):

    c = lambda s: float(s.decode().replace('D', 'e'))
    N = 9
    data = np.genfromtxt(file_name, delimiter=",", converters=dict(zip(range(N), [c]*N)))

    for ID in np.arange(8):
        plt.figure(ID)
        plt.plot(data[:, 0], data[:, ID+1], label=legend)
        plt.xlabel('t [s]')
        plt.grid(True)
        title_fig = title + " - fracture id " + str(ID)
        plt.title(title_fig)
        plt.legend()
        plt.ylabel(r"$\\\overline{c_2}$")

def save(filename):
    folder = "./plots/"

    if not os.path.exists(folder):
        os.makedirs(folder)

    for ID in np.arange(8):
        plt.figure(ID)
        plt.savefig(folder+filename+"_fracture_"+str(ID), bbox_inches='tight')
        plt.gcf().clear()

#------------------------------------------------------------------------------#
# add data to plots

# add plots for the different schemes
# The first argument is the path to your data, where the data is assumed to be ordered comma-separated in the following order:
#   -> 1. time
#      2. integral of \phi c within matrix sub-domain \Omega_3
#      3. integral of \phi c within fracture domain \Omega_f
#      4. outlux of the domain across the outfow boundary
# The second argument defines the legend you want to add to your data (Institution / Numerical method)
# The third argument specifies the plot id - use the ids defined in lines 35-37 for the different plots

# TODO: add reference solution to plots as soon as available

titles = ['$\\sim 30k$ cells', '$\\sim 150k$ cells']
refinement_index = ['0', '1']

for title, ref in zip(titles, refinement_index):

#    # University of Stuttgart Mpfa results
#    data = "../results/USTUTT/MPFA/dataovertime_"+ref+".txt"
#    plot(data, "USTUTT-MPFA", id_intc_matrix, title)
#    plot(data, "USTUTT-MPFA", id_intc_fracture, title)
#    plot(data, "USTUTT-MPFA", id_outflux, title)

    # University of Bergen results
    methods = ['MVEM', 'TPFA', 'MPFA', 'RT0']
    for method in methods:
        data = "../results/UiB/"+method+"/dot_refinement_"+ref+".csv"
        label = "UiB-"+method
        plot(data, label, title)

    # UNICE-UNIGE
    methods = ['HFVCont', 'HFVDisc', 'VAGCont', 'VAGDisc']
    for method in methods:
        data = "../results/UNICE-UNIGE/"+method+"/dot_refinement_"+ref+".csv"
        label = "UNICE-UNIGE-"+method
        plot(data, label, title)

    # save figures
    save("pot_refinement_"+ref)

#------------------------------------------------------------------------------#
