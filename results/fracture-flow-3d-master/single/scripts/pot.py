import matplotlib.pyplot as plt
import numpy as np
import os
import sys

#------------------------------------------------------------------------------#

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=15)

# ids of the different plots
id_intc_matrix = 0   # integral of c*porosity in the lower matrix sub-domain
id_intc_fracture = 1 # indegral of c*porosity in the fracture
id_outflux = 2       # integrated outflux across the outflow boundary

def plot(file_name, legend, ID, title):

    c = lambda s: float(s.decode().replace('D', 'e'))
    N = 4
    data = np.genfromtxt(file_name, delimiter=",", converters=dict(zip(range(N), [c]*N)))

    plt.figure(ID)
    plt.plot(data[:, 0]/(365*24*3600), data[:, ID+1], label=legend)
    plt.xlabel('t [y]')
    plt.grid(True)
    plt.title(title)
    if ID == id_intc_fracture:
        plt.legend(loc="lower right")
    else:
        plt.legend(loc="upper left")

    # choose y-label depending on plot id
    if ID == id_intc_matrix:
        plt.ylabel(r"$\int_{\Omega_3} \phi c$")
    elif ID == id_intc_fracture:
        plt.ylabel(r"$\int_{\Omega_f} \phi c$")
    elif ID == id_outflux:
        plt.ylabel("outflux")
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
#   -> 1. time
#      2. integral of \phi c within matrix sub-domain \Omega_3
#      3. integral of \phi c within fracture domain \Omega_f
#      4. outlux of the domain across the outfow boundary
# The second argument defines the legend you want to add to your data (Institution / Numerical method)
# The third argument specifies the plot id - use the ids defined in lines 35-37 for the different plots

# TODO: add reference solution to plots as soon as available

titles = ['$\\sim 1k$ cells', '$\\sim 10k$ cells', '$\\sim 100k$ cells']
refinement_index = ['0', '1', '2']

for title, ref in zip(titles, refinement_index):

    # University of Stuttgart Mpfa results
    data = "../results/USTUTT/MPFA/dot_refinement_"+ref+".txt"
    plot(data, "USTUTT-MPFA", id_intc_matrix, title)
    plot(data, "USTUTT-MPFA", id_intc_fracture, title)
    plot(data, "USTUTT-MPFA", id_outflux, title)

    # University of Bergen results
    methods = ['MVEM', 'TPFA', 'MPFA', 'RT0']
    for method in methods:
        data = "../results/UiB/"+method+"/dot_refinement_"+ref+".csv"
        label = "UiB-"+method
        plot(data, label, id_intc_matrix, title)
        plot(data, label, id_intc_fracture, title)
        plot(data, label, id_outflux, title)

    # UNICE-UNIGE
    methods = ['HFVCont', 'HFVDisc', 'VAGCont', 'VAGDisc']
    for method in methods:
        data = "../results/UNICE-UNIGE/"+method+"/dot_refinement_"+ref+".csv"
        label = "UNICE-UNIGE-"+method
        plot(data, label, id_intc_matrix, title)
        plot(data, label, id_intc_fracture, title)
        plot(data, label, id_outflux, title)

    # NCU_TW
    methods = ["NAME"]
    for method in methods:
        data = "../results/NCU_TW/"#+method
        data += "/dot_refinement_"+ref+".csv"
        label = "NCU-TW-"+method
        plot(data, label, id_intc_matrix, title)
        plot(data, label, id_intc_fracture, title)
        plot(data, label, id_outflux, title)

    # save figures
    save(id_intc_matrix, "pot_c_matrix"+ref)
    save(id_intc_fracture, "pot_c_fracture"+ref)
    save(id_outflux, "pot_outflux"+ref)

#------------------------------------------------------------------------------#
