import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import matplotlib.animation as animation

plt.rcParams["font.family"] = 'Arial'
# plt.rcParams["font.sans-serif"] = "DejaVu Sans, Arial, Helvetica, Lucida Grande, Verdana, Geneva, Lucid, Avant Garde, sans-serif"
plt.rcParams["font.size"] = 12
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["xtick.major.size"] = 4.75
plt.rcParams["xtick.major.width"] = 0.75
plt.rcParams["xtick.minor.size"] = 3
plt.rcParams["xtick.minor.width"] = 0.4
plt.rcParams["xtick.minor.visible"] = True

plt.rcParams["ytick.direction"] = "in"
plt.rcParams["ytick.major.size"] = 4.75
plt.rcParams["ytick.major.width"] = 0.75
plt.rcParams["ytick.minor.size"] = 3
plt.rcParams["ytick.minor.width"] = 0.4
plt.rcParams["ytick.minor.visible"] = True

plt.rcParams["savefig.transparent"] = True

plt.rcParams['axes.linewidth'] = 0.75

rootdir = "./"
# rootdir = "../debug_build/example_4mNonSelfSim"

simulation_id = "fb03-087__0129_a=4.00_ruptype=pulse_pdcscaling=0.600_sn=6.0MPa_hatfr=0.3_bgbeta=0.35"
field_name = "top_velo_0" # select plot parameter

# read meta data of simulation
df_time = pd.read_csv(os.path.join(rootdir, simulation_id+".time"), header=None, sep=' ', index_col=0)
df_coord = pd.read_csv(os.path.join(rootdir, simulation_id+".coords"), header=None, sep=' ', index_col=None)
NT=len(df_time)

xcoord = df_coord.loc[:,0].values
zcoord = df_coord.loc[:,2].values

x_length = xcoord.max()
z_length = zcoord.max()

nb_x_elements = int(np.sqrt(len(df_coord)))
nb_z_elements = nb_x_elements # assume same number of elements
# ref also: https://stackoverflow.com/a/35176314

X = xcoord.reshape(nb_x_elements,nb_z_elements).T - x_length/2
Z = zcoord.reshape(nb_x_elements,nb_z_elements).T - z_length/2


# read data
IfBinaryOutput = True # true if the output is in binary format
if IfBinaryOutput:
    D = np.fromfile(os.path.join(rootdir,simulation_id+f"-DataFiles/{field_name}.out"), dtype="float32")
    df_data = pd.DataFrame(data=D.reshape((NT, -1)))
else:
    df_data = pd.read_csv(os.path.join(rootdir,simulation_id+f"-DataFiles/{field_name}.out"), header=None, sep=' ', engine="c")
    
if "top_" in field_name:
    value_factor_double = 2.0 # double for slip and slip velocity
else:
    value_factor_double = 1.0

# Plot animation
fig, ax = plt.subplots(1, 1, figsize=(7, 6))

ims=[]

plotstep=1 # decimate the snap shots when setting plotstep > 1

vmin, vmax = 0, np.max(df_data)

for i, tplot in enumerate(df_time[1].values[1::plotstep]):
    
    t_plot_ind = np.where(df_time[1].values >= tplot)[0][0]
    data_snap = df_data.loc[t_plot_ind, :].values

    # print(t_plot_ind, tplot)
    
    V = data_snap.reshape(nb_x_elements,nb_z_elements).T    
    h1 = ax.pcolormesh(X*1e3, Z*1e3, value_factor_double*V, cmap='inferno', rasterized=True, vmin=vmin, vmax=vmax) # multiply 2 for slip velocity
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("z [mm]")
    ax.set_aspect('equal', 'box')
    # fig.colorbar(h1, orientation="vertical", shrink=1.0)

    ax.set_xlim([-12, 12])
    ax.set_ylim([-12, 12])

    # ax.set_title does not work for the animation
    titlestr = f"{field_name} T={tplot*1e6:.2f}us"
    title = ax.text(0.5, 1.01, titlestr, horizontalalignment='center', verticalalignment='bottom', transform=ax.transAxes)

    if i == 0:
        fig.colorbar(h1, ax=ax, orientation="vertical",  shrink=0.9, label='Velocity [m/s]')
            
    ims.append([h1, title])

ani = animation.ArtistAnimation(fig, ims, blit=False, interval=15, repeat_delay=2000)

# Save animation
ani.save(f'./animation_{field_name}.gif', writer='ffmpeg', fps=30, dpi=100)

plt.show()

