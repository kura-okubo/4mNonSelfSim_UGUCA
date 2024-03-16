#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

exname = 'fracture_2d_example_restarted'

def plot(fldname,comp,ax):

    # get time data
    Tdata = []
    with open(exname+".time",'r') as fl:
        lines = fl.readlines()
    for line in lines:
        Tdata.append(float(line.strip().split()[-1]))

    # get space data
    Xdata = []
    with open(exname+".coords",'r') as fl:
        lines = fl.readlines()
    for line in lines:
        Xdata.append(float(line.strip().split()[0]))

    # get field data
    Vdata = []
    with open(exname+"-DataFiles/"+fldname+".out",'r') as fl:
        lines = fl.readlines()
    for line in lines:
        Vdata.append([float(i) for i in line.strip().split()])
    Vdata = np.array(Vdata)

    # get components
    nb_nodes = (int)(Vdata.shape[1]/2)
    Vdata = Vdata[:,comp*nb_nodes:(comp+1)*nb_nodes]

    # plot
    XV, TV = np.meshgrid(Xdata, Tdata)
    pc = ax.pcolor(XV,TV,Vdata)
    return pc

# -----------------------------------------------------------------------------
if __name__ == "__main__":

    if len(sys.argv) not in [3]:
        sys.exit('Missing argument! usage: ./fracture_2d_example_plot.restart.py '
                 + 'fieldname '
                 + '(options: cohesion top_disp bot_disp tau_max)'
                 + 'comp (options: 0 1)')

    fldname=str(sys.argv[1])
    comp=int(sys.argv[2])
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pc = plot(fldname,comp,ax)
    cbar = fig.colorbar(pc)
    cbar.set_label(fldname)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    plt.show()
