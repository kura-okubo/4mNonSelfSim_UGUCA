#!/usr/bin/env python3
from __future__ import print_function
from __future__ import division
from glob import glob
import os
import sys
import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from plotting_utils import *



def generate_output(full_path, bname):
  stations = ['fltst_dp000', 'fltst_dp025', 'fltst_dp050', 
              'fltst_dp075', 'fltst_dp100', 'fltst_dp125',  
              'fltst_dp150', 'fltst_dp175', 'fltst_dp200',  
              'fltst_dp250', 'fltst_dp300', 'fltst_dp350']
  if not os.path.exists('upload'):
    os.makedirs('upload')
  for station in stations:
    dump_station(full_path, bname, station)
  # dump_cplot(full_path, bname)


# def dump_cplot(full_path, bname):
#   dump_path = 'upload/%s' % bname
#   if not os.path.exists(dump_path):
#     os.makedirs(dump_path)
#   nb_nodes_x = int(bname.split('_Nx')[1].split('_')[0])
#   domain_factor = float(bname.split('_s')[1].split('_')[0])
#   nb_nodes = nb_nodes_x * nb_nodes_z
#   length_x_rpt = 36e3
#   length_z_rpt = 18e3
#   length_x = domain_factor * length_x_rpt
#   length_z = domain_factor * length_z_rpt
#   dx = length_x / nb_nodes_x
#   dz = length_z / nb_nodes_z
#   x = np.arange(nb_nodes_x) * dx - length_x / 2
#   z = 7.5e3 + length_z / 2 - np.arange(nb_nodes_z) * dz
#   t = np.fromfile('%s.time' % full_path, sep=' ')
#   t = t[1:-1:2]
#   nt = len(t)

#   delta_dot = read_data_cplot('%s-DataFiles/top_velo_0.out' % full_path, nb_nodes_x, nb_nodes_z, nt) * 2

#   C = np.zeros((nb_nodes_x * nb_nodes_z, 3))
#   for j in range(nb_nodes_x):
#     for k in range(nb_nodes_z):
#       i = np.argmax(delta_dot[:, j, k] > 1e-3)
#       C[j * nb_nodes_z + k, 0] = x[j]
#       C[j * nb_nodes_z + k, 1] = z[k]
#       C[j * nb_nodes_z + k, 2] = t[i]
#   idx = np.argsort(C[:, 2], kind='mergesort')
#   C = C[idx]

#   with open('%s/ke_%s.txt' % (dump_path, 'cplot'), 'w') as f:
#     meta = [
#       '# problem = TPV101',
#       '# author = Kammer, Albertini, Ke',
#       '# date = 2021/08/20',
#       '# code = uguca',
#       '# code_version = 0.9',
#       '# element_size = %d m' % dx,
#       '# time_step = %.5e s' % (t[1] - t[0]),
#       '# num_time_steps = %d' % nt,
#       '# Column #1 = horizontal coordinate, distance along strike (m)',
#       '# Column #2 = vertical coordinate, distance down-dip (m)',
#       '# Column #3 = rupture time (s)']
#     f.writelines('\n'.join(meta))
#     f.write('# \n')
#     f.write('j  k  t\n')
#     for i in range(C.shape[0]):
#       if C[i, 0] < -15000 or C[i, 0] > 15000 or C[i, 1] > 15000 or C[i, 1] < 0:
#         continue
#       f.write('%14.6e%14.6e%14.6e\n' % tuple(C[i, :]))


def dump_station(full_path, bname, station):
  dump_path = 'upload/%s' % bname
  if not os.path.exists(dump_path):
    os.makedirs(dump_path)
  location = station.split('fltst')[1].split('dp')
  x_interest = int(location[1]) * 100
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
  print('Warning: %.2e = %.2e? ' % (x_interest, x[idx_x]))

  t = np.fromfile('%s.time' % full_path, sep=' ')
  t = t[1::2]
  nt = len(t)
  delta = read_data('%s-DataFiles/top_disp.out' % full_path, nb_nodes_x, 1, nt, idx_x, 0) * 2
  delta = delta[:, 0]
  delta_dot = read_data('%s-DataFiles/top_velo.out' % full_path, nb_nodes_x, 1, nt, idx_x, 0) * 2
  delta_dot = delta_dot[:, 0]
  cohesion = read_data('%s-DataFiles/cohesion.out' % full_path, nb_nodes_x, 1, nt, idx_x, 0)
  cohesion = cohesion[:, 0]
  theta = read_data('%s-DataFiles/theta.out' % full_path,
                    nb_nodes_x, 1, nt, idx_x, 0)
  with open('%s/ke_%s.txt' % (dump_path, station), 'w') as f:
    meta = [
      '# problem = BP1-FD',
      '# author = Kammer, Albertini, Ke',
      '# date = 2024/03/17',
      '# code = uguca',
      '# code_version = 1.0', 
      '# element_size = %d m' % dx,
      '# time_step = %.5e s' % (t[1] - t[0]),
      '# num_time_steps = %d' % nt,
      '# location= on fault %.2f km along strike, %.2f km down-dip' % (x_interest / 1000, z_interest / 1000),
      '# Column #1 = Time (s)',
      '# Column #2 = horizontal slip (m)',
      '# Column #3 = log horizontal slip rate (m/s)',
      '# Column #4 = horizontal shear stress (MPa)',
      '# Column #4 = log theta']
    f.writelines('\n'.join(meta))
    f.write('# \n')
    f.write('t  h-slip  h-slip-rate  h-shear-stress  log-theta\n')
    
    for i in range(nt):
      f.write(('%20.12e' + '%14.6e' * 5 + '\n') % (t[i], delta[i], np.log10(delta_dot[i]), cohesion[i] * 1e-6, np.log10(theta[i])))


if __name__ == '__main__':
  usage = """./inspect_results.py <path_to_bname>"""
  if len(sys.argv) != 2:
      sys.exit(usage)
  try:
    full_path = sys.argv[1]
    bname = sys.argv[1].split('/')[-1]
  except:
    sys.exit(usage)

  generate_output(full_path, bname)
