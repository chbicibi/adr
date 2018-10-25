# -*- coding: utf-8 -*-
import sys, os, glob, re, subprocess
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import csv

def plot_ball(ax):
  r = 6378
  u = np.linspace(0, 2 * np.pi, 100)
  v = np.linspace(0, np.pi, 100)
  x = r * np.outer(np.cos(u), np.sin(v))
  y = r * np.outer(np.sin(u), np.sin(v))
  z = r * np.outer(np.ones(np.size(u)), np.cos(v))
  ax.plot_wireframe(x, y, z, rstride=6, cstride=6, color='gray', alpha=0.2)
  return x, y, z

def plot_orbit(ax, trj, label):
  # r  = 3
  # nu = np.linspace(0, 2 * np.pi, num=30)
  # x  = r * np.cos(nu)
  # y  = r * np.sin(nu)
  # z  = 1 * np.cos(nu + 0.5 * np.pi)
  x, y, z = trj
  ax.plot(x, y, z, label=label)
  return x, y, z

def load_orbit(file):
  if not os.path.exists(file):
    print("入力ファイルがありません: " + file)
    return

  # labels = ["$\\alpha, deg$", "$elv, deg$", "$Ma$", "$x_1$", "$x_2$"]
  # x_label = labels[3]
  # y_label = labels[4]
  # z_label = ["$c_x$", "$c_m$", "$c_z$", "$z$"][int(no) - 1]

  with open(file, "r") as f:
    x, y, z = (np.array(x, dtype=np.float) for x in list(zip(*(x for x in csv.reader(f) if x))))
  return (x, y, z)

def transparent(ifile, ofile):
  cmd = "convert %s -transparent white %s" % (ifile, ofile)
  print(cmd)
  returncode = subprocess.call(cmd)

def main(argv):
  fig = plt.figure(figsize=(8, 7))
  fig.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.0, wspace=0.0, hspace=0.0)
  # ax  = Axes3D(fig)
  ax  = fig.add_subplot(111, projection='3d', aspect='equal')
  plot_ball(ax)

  trj = load_orbit("orbit1.csv")
  plot_orbit(ax, trj, label="departure")

  trj = load_orbit("orbit2.csv")
  plot_orbit(ax, trj, label="arrival")

  trj = load_orbit("orbit3.csv")
  plot_orbit(ax, trj, label="transfer")

  plt.axis('off')
  plt.legend(loc='best', fontsize=14, frameon=False)

  lim = 6000
  ax.set_xlim(-lim, lim)
  ax.set_ylim(-lim, lim)
  ax.set_zlim(-lim, lim)

  ax.view_init(30, 240)
  # ax.view_init(72.5, 240)
  # ax.set_aspect([1, 1, 1])

  if len(argv) == 0:
    plt.show()
  else:
    file = argv[0]
    fig.savefig(file + ".png", dpi=100, transparent=True)

main(sys.argv[1:])
# transparent("Figure_6.png", "Figure_6t.png")
