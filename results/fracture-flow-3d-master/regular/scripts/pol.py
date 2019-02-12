import matplotlib.pyplot as plt
import numpy as np
import os

#------------------------------------------------------------------------------#

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=15)

def plot(file_name, ID, simulation_id, title, cond):

    c = lambda s: float(s.decode().replace('D', 'e'))
    N = 2
    data = np.genfromtxt(file_name, delimiter=",", converters=dict(zip(range(N), [c]*N)))

    plt.figure(simulation_id)
    plt.plot(data[:, 0], data[:, 1], label=ID)
    plt.xlabel('archlenght')
    plt.ylabel('$p$')
    plt.grid(True)
    plt.xlim((0, np.sqrt(3)))
    if cond == 0:
        plt.ylim((0.95, 3))
        title += " - permeability 1e4"
    else:
        plt.ylim((0, 15.5))
        title += " - permeability 1e-4"
    plt.legend()
    plt.title(title)

def save(simulation_id, filename):
    folder = "./plots/"

    if not os.path.exists(folder):
        os.makedirs(folder)

    plt.figure(simulation_id)
    plt.savefig(folder+filename, bbox_inches='tight')
    plt.gcf().clear()

#------------------------------------------------------------------------------#
# Insert here your data for the plotting, see the file 'color_regions.vtu'
# for the coloring code of each region.

titles = ['$\\sim 500$ cells', '$\\sim 4k$ cells', '$\\sim 32k$ cells']
refinement_index = [0, 1, 2]
conds = [0, 1]

for cond in conds:
    for title, ref in zip(titles, refinement_index):

        # University of Bergen
        methods = ['MVEM', 'TPFA', 'MPFA', 'RT0']
        for method in methods:
            folder = "../results/UiB/"+method+"/"
            label = "UiB-"+method
            name = folder+"dol_cond_"+str(cond)+"_refinement_"+str(ref)+".csv"
            plot(name, label, ref, title, cond)

        # University of Stuttgart
#        folder = "../results/USTUTT/MPFA/"
#        label = "USTUTT-MPFA"
#        name = folder+"dol_perm_"+str(cond)+"_refinement_"+str(ref)+".csv"
#        plot(name, label, ref, title, cond)

        # UNICE-UNIGE
        methods = ['HFVCont', 'HFVDisc', 'VAGCont', 'VAGDisc']
        for method in methods:
            folder = "../results/UNICE-UNIGE/"+method+"/"
            label = "UNICE-UNIGE-"+method
            name = folder+"dol_cond_"+str(cond)+"_refinement_"+str(ref)+".csv"
            plot(name, label, ref, title, cond)

        save(ref, "pol_cond_"+str(cond)+"_refinement_"+str(ref))

#------------------------------------------------------------------------------#
