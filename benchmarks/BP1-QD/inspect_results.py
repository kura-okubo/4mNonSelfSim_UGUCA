#!/usr/bin/env python3
from __future__ import print_function
from __future__ import division
from glob import glob
import sys
import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from plotting_utils import *

dont_exclude_author = np.array(['dalzillo', 'jiang', 'liang', 'yang'])

dont_exclude_author_c = np.array(['dalzillo', 'jiang', 'liang', 'yang'])


def inspect_results(full_path, bname):
  stations = glob('./ref/*dp*/')
  if not stations:
    stations = ['fltst_dp000/', 'fltst_dp075/',
                'fltst_dp150/', 'fltst_dp200/']
  for station in stations:
    location = station.split('fltst')[1].split('/')[-2]
    location = location.split('dp')
    x = int(location[1]) * 100
    compare_station(full_path, bname, x)
  # compare_cplot(full_path, bname)


# def compare_cplot(full_path, bname):
#   nb_nodes_x = int(bname.split('_Nx')[1].split('_')[0])
#   nb_nodes_z = int(bname.split('_Nz')[1].split('_')[0])
#   domain_factor = float(bname.split('_s')[1].split('_')[0])
#   nb_nodes = nb_nodes_x * nb_nodes_z
#   length_x_rpt = 36e3
#   length_z_rpt = 18e3
#   length_x = domain_factor * length_x_rpt
#   length_z = domain_factor * length_z_rpt
#   dx = length_x / nb_nodes_x
#   dz = length_z / nb_nodes_z
#   x = np.arange(nb_nodes_x) * dx - length_x / 2
#   z = -7.5e3 + length_z / 2 - np.arange(nb_nodes_z) * dz
#   t = np.fromfile('%s.time' % full_path, sep=' ')
#   t = t[1:-1:2]
#   nt = len(t)
#   xx, zz = np.meshgrid(x, z, indexing='ij')

#   delta_dot = read_data_cplot('%s-DataFiles/top_velo_0.out' %
#                          full_path, nb_nodes_x, nb_nodes_z, nt) * 2

#   rpt_arrival = np.zeros(xx.shape)
#   for j in range(nb_nodes_x):
#     for k in range(nb_nodes_z):
#       i = np.argmax(delta_dot[:, j, k] > 1e-3)
#       rpt_arrival[j, k] = t[i]

#   fig, ax = plt.subplots(1)
#   times = np.arange(0, 8, 0.5)

#   stname = 'cplot'
#   data_dir = './ref/'+stname+'/*'
#   refs = glob(data_dir)
#   i = 0
#   for ref in refs:
#     author = ref.split('/')[-1].split('.txt')[0]
#     if author in dont_exclude_author_c:
#       c = 'bgrcmy'[i]
#       i += 1
#       pass
#     else:
#       continue
#     X, Y, T = read_scec_cplot(ref)  # ,nb_nodes_x,nb_nodes_z,nt)
#     ax.contour(X, Y, T, colors=c, linestyles='--', levels=times)
#     ax.plot([], [], color=c, linestyle='--', label=author)
#   ax.contour(xx, zz, rpt_arrival, colors='k', levels=times)
#   ax.plot([], [], color='k', linestyle='-', label='uguca')
#   #ax.pcolor(xx.T,zz.T,tau0[200])

#   ax.legend(loc='upper right', ncol=1, fontsize=7)

#   ax.set_ylim(-7.5e3 - length_z_rpt / 2, -7.5e3 + length_z_rpt / 2)
#   ax.set_xlim(-length_x_rpt/2, length_x_rpt/2)
#   ax.set_aspect(1)
#   ax.set_xlabel(r'$x$ (m)')
#   ax.set_ylabel(r'$y$ (m)')
#   plt.tight_layout()
#   print('savefig {}_{}.pdf'.format(bname, stname))
#   plt.savefig('{}_{}.pdf'.format(bname, stname))


def compare_station(full_path, bname, x_interest):
  spec = bname.split('_')
  nb_nodes_x = int(spec[1][1::])
  domain_factor = float(spec[2][1::])
  nb_nodes = nb_nodes_x
  length_x_rpt = 50e3
  length_x = domain_factor * length_x_rpt
  dx = length_x / nb_nodes_x
  x = np.arange(nb_nodes_x) * dx - length_x / 2
  idx_x = np.argmin(np.abs(x - x_interest))
  idx = idx_x
  print('Double check: %.2e = %e ?' % (x_interest, x[idx_x]))

  t = np.fromfile('%s.time' % full_path, sep=' ')
  t = t[1:-10:2]
  nt = len(t)
  delta = read_data('%s-DataFiles/top_disp.out' % full_path, nb_nodes_x, 1, nt, idx_x, 0) * 2
  delta = delta[:, 0]
  delta_dot = read_data('%s-DataFiles/top_velo.out' % full_path, nb_nodes_x, 1, nt, idx_x, 0) * 2
  delta_dot = delta_dot[:, 0]
  cohesion = read_data('%s-DataFiles/cohesion.out' % full_path, nb_nodes_x, 1, nt, idx_x, 0)
  cohesion = cohesion[:, 0]
  theta = read_data('%s-DataFiles/theta.out' % full_path, nb_nodes_x, 1, nt, idx_x, 0)

  params = {
      'text.latex.preamble': r'''\usepackage{siunitx}',
                                 \usepackage{sfmath}
                                 \sisetup{detect-family = true}
                                 \usepackage{amsmath}''',
      'font.family': 'sans-serif',
      'font.sans-serif': ['Helvetica']
  }
  plt.rcParams.update(params)
  fig = plt.figure()
  fig.subplots_adjust(top=0.97, bottom=0.1, left=0.1,
                      right=0.97, wspace=0.27, hspace=0.35)
  ax1 = fig.add_subplot(2, 2, 1)
  ax1.plot(t, delta, '-k', label='uguca', linewidth=1, zorder=10)
  ax1.set_xlabel(r'$t$ (sec)')
  ax1.set_ylabel(r'$\delta$ (m)')
  ax2 = fig.add_subplot(2, 2, 2)
  ax2.semilogy(t, delta_dot, '-k', label='uguca', linewidth=1, zorder=10)
  ax2.set_xlabel(r'$t$ (sec)')
  ax2.set_ylabel(r'$\dot{\delta}$ (m/s)')
  ax3 = fig.add_subplot(2, 2, 3)
  ax3.plot(t, cohesion * 1e-6, '-k', label='uguca', linewidth=1, zorder=10)
  ax3.set_xlabel(r'$t$ (sec)')
  ax3.set_ylabel(r'$\tau$ (MPa)')
  ax4 = fig.add_subplot(2, 2, 4)
  ax4.semilogy(t, theta, '-k', label='uguca', linewidth=1, zorder=10)
  ax4.set_xlabel(r'$t$ (sec)')
  ax4.set_ylabel(r'$\theta$ (m/s)')

  data_dir = './ref/fltst_dp%03d/*' % (x_interest / 100)
  refs = glob(data_dir)
  for ref in refs:
    author = ref.split('/')[-1].split('.')[0]
    data = read_scec(ref)
    ax1.plot(data['time'], data['delta_x'], '--', label=author, linewidth=1)
    ax2.semilogy(data['time'], 10 ** data['delta_dot_x'],
             '--', label=author, linewidth=1)
    ax3.plot(data['time'], data['tau_x'], '--', label=author, linewidth=1)
    ax4.semilogy(data['time'], 10 ** data['log_theta'],
                 '--', label=author, linewidth=1)

  ax2.legend(loc='upper right', ncol=1, fontsize=7)
  plt.savefig('%s_dp%03d.pdf' %(bname, x_interest / 100))


def read_scec(path):
  raw = np.genfromtxt(path, skip_header=1)
  raw = np.delete(raw, 0, 0)  # don't know why it didn't discard column headers
  raw = np.array(raw)
  data = dict()
  data['time'] = raw[1::, 0]
  data['delta_x'] = raw[1::, 1]
  data['delta_dot_x'] = raw[1::, 2]
  data['tau_x'] = raw[1::, 3]
  data['log_theta'] = raw[1::, 4]
  return data


def read_scec_cplot(path):
  raw = np.genfromtxt(path, skip_header=1)
  raw = np.delete(raw, 0, 0)  # don't know why it didn't discard column headers
  raw = np.array(raw)
  X = raw[0::, 0]  # distance along strike
  Y = raw[0::, 1]  # distance along dip
  Z = raw[0::, 2]  # rpture time (s)
  ny = len(X[X == X[0]])
  nx = len(X[Y == Y[0]])
  n = len(X)
  if n != ny*nx:
    raise RuntimeError("nx*ny!=n")
  if True:
    if X[0] == X[1]:
      X = X.reshape(nx, ny)
      Y = -Y.reshape(nx, ny)
      Z = Z.reshape(nx, ny)
    else:
      x = np.linspace(np.min(X), np.max(X), 100)
      y = np.linspace(np.min(Y), np.max(Y), 100)
      triang = tri.Triangulation(X, Y)
      interpolator = tri.LinearTriInterpolator(triang, Z)
      Xi, Yi = np.meshgrid(x, y)
      Zi = interpolator(Xi, Yi)
      X = Xi
      Y = -Yi
      Z = Zi
    return X, Y, Z


if __name__ == '__main__':
  usage = """./inspect_results.py <path_to_bname>"""
  if len(sys.argv) != 2:
      sys.exit(usage)
  try:
    full_path = sys.argv[1]
    bname = sys.argv[1].split('/')[-1]
  except:
    sys.exit(usage)

  inspect_results(full_path, bname)
