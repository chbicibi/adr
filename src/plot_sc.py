# coding: utf-8
import sys, os, glob, re, subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
import csv


def plotf2d(X, Y, Z, lx, ly, lz, file=None):
  fig = plt.figure(figsize=(8, 7))
  ax = plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.95, wspace=0.2, hspace=0.2)
  # ax = Axes3D(fig)
  plt.xlabel(lx, size=14)
  plt.ylabel(ly, size=14)
  plt.tick_params(labelsize=14)

  # im = plt.pcolor(X, Y, Z, cmap=None, linewidth=0, vmin=0, vmax=50)
  plt.scatter(X, Y, s=10, c=Z, cmap="Blues", vmin=-30, vmax=100)

  plt.xlim(0, 50)
  plt.ylim(0, 0.4)
  # ax.set_zlim(0, 100)
  # cbar = fig.colorbar(im)
  # cbar.set_label(lz, size=14)
  # plt.contour(X, Y, Z, zdir='z')
  if file:
    plt.savefig(file)
  else:
    plt.show()
  # ax.imwrite("out.png")

def plotf3d(X, Y, Z, lx, ly, lz):
  fig = plt.figure()
  fig.subplots_adjust(left=0.15, bottom=0.1, right=0.7, top=0.95, wspace=0.1, hspace=0.2)
  ax = Axes3D(fig)
  ax.set_xlabel(lx, size=14)
  ax.set_ylabel(ly, size=14)
  ax.set_zlabel(lz, labelpad=30, size=14)
  ax.tick_params(labelsize=14)
  ax.tick_params(axis='z', pad=20)

  # plt.gca().zaxis.set_tick_params(which='both', direction='in',bottom=True, top=True, left=True, right=True)
  # ax.set_aspect(0.2)

  ax.plot_surface(X, Y, Z, cmap='bwr', linewidth=0)
  ax.contour(X, Y, Z, zdir='z', offset=np.min(Z))
  plt.show()
  # ax.imwrite("out.png")

def plotf3d2(X, Y, Z, P, Q, R, lx, ly, lz):
  fig = plt.figure()
  fig.subplots_adjust(left=0.15, bottom=0.1, right=0.7, top=0.95, wspace=0.1, hspace=0.2)
  ax = Axes3D(fig)
  ax.set_xlabel(lx, size=14)
  ax.set_ylabel(ly, size=14)
  ax.set_zlabel(lz, labelpad=30, size=14)
  ax.tick_params(labelsize=14)
  ax.tick_params(axis='z', pad=20)

  # plt.gca().zaxis.set_tick_params(which='both', direction='in',bottom=True, top=True, left=True, right=True)
  # ax.set_aspect(0.2)

  ax.plot_surface(X, Y, Z, cmap='bwr', linewidth=0)
  ax.plot(P, Q, R, "o")
  ax.contour(X, Y, Z, zdir='z', offset=np.min(Z))
  plt.show()
  # ax.imwrite("out.png")

def main(argv):
  files = glob.glob('res*.csv')
  if len(argv) == 0:
    # print("種類を指定してください: nsga2=>0, moead=>1")
    # ty = sys.stdin.readline().rstrip("\n")
    print("ファイル番号を指定してください:")
    print("[-1] なし")
    for i, f in enumerate(files):
      print(f"[{i}] {f}")
    no = int(sys.stdin.readline())
  elif len(argv) == 1:
    no = int(argv[0])
  else:
    return

  if not no in range(len(files)):
    return
  file = files[no]

  if not os.path.exists(file):
    print(f"入力ファイルがありません: {file}")
    return
  else:
    print(f"入力ファイル: {file}")

  labels = ["$\\Delta V, km/s$", "$rcs, m^2$", "$f$"]
  x_label = labels[0]
  y_label = labels[1]
  z_label = labels[2]

  ofile = "out_nsga2.png"

  with open(file, "r") as f:
    reader = csv.reader(f)
    next(reader)
    Z, X, Y = (np.array(x, dtype=np.float) for x in list(zip(*reader))[0:3])
  plotf2d(X, Y, Z, x_label, y_label, z_label, file=ofile)

main(sys.argv[1:])
