import matplotlib.pyplot as plt
import numpy as np
import os

#------------------------------------------------------------------------------#

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=15)

main_folder = "./plots/"

def plot(file_name, legend, title, num_colors = 21):

    c = lambda s: float(s.decode().replace('D', 'e'))
    N = 22
    data = np.genfromtxt(file_name, delimiter=",", converters=dict(zip(range(N), [c]*N)))

    for color in np.arange(num_colors):
        plt.figure(color)
        plt.plot(data[:, 0], data[:, color+1], label=legend)
        plt.title(title)
        plt.xlabel('$t$')
        plt.ylabel('$c$')
        plt.grid(True)
        plt.legend()

#------------------------------------------------------------------------------#

def save(filename, cond, num_colors = 21):

    if not os.path.exists(main_folder):
        os.makedirs(main_folder)

    folder = main_folder + "cond_" + str(cond) + "/"

    if not os.path.exists(folder):
        os.makedirs(folder)

    for color in np.arange(num_colors):
        plt.figure(color)
        name = filename + "_region_" + str(color)
        plt.savefig(folder+name, bbox_inches='tight')
        plt.gcf().clear()

#------------------------------------------------------------------------------#

# Insert here your data for the plotting, see the file 'color_regions.vtu'
# for the coloring code of each region.


titles = ['$\\sim 4k$ cells  - permeability 1e4', '$\\sim 4k$ cells  - permeability 1e-4']
conds = [0, 1]

#------------------------------------------------------------------------------#

for cond, title in zip(conds, titles):
    # University of Bergen results
    methods = ['MVEM', 'TPFA', 'MPFA', 'RT0']
    for method in methods:
        data = "../results/UiB/"+method+"/dot_cond_"+str(cond)+".csv"
        label = "UiB-"+method
        plot(data, label, title)

#    # university of stuttgart results
#    data = "../results/USTUTT/MPFA/dot_perm_"+str(cond)+".csv"
#    label = "USTUTT-MPFA"
#    plot(data, label, title)

    # UNICE-UNIGE
    methods = ['HFVCont', 'HFVDisc', 'VAGCont', 'VAGDisc']
    for method in methods:
        data = "../results/UNICE-UNIGE/"+method+"/dot_cond_"+str(cond)+".csv"
        label = "UNICE-UNIGE-"+method
        plot(data, label, title)

    name = "cot_cond_"+str(cond)
    save(name, cond)

#------------------------------------------------------------------------------#
