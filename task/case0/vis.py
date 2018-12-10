#! /usr/bin/env python3

'''
Abstruct
'''

import argparse
import os
import shutil

import numpy as np
import matplotlib.pyplot as plt

import orbit_model as case


################################################################################


def plot_hv(out='result', show=False):
    dir_res = 'result'
    files = [
        ('NSGA-II', 'nsga2/pop_size500_epoch100_20181121_1010_hv.npy'),
        ('MOEA/D case1', 'moead/pop_size500_epoch100_20181121_2320_hv.npy'), # all
        ('MOEA/D case2', 'moead/pop_size500_epoch100_20181121_1006_hv.npy'), # new
        # ('MOEA/D new', 'moead/pop_size500_epoch100_20181122_0425_hv.npy'),
    ]

    fig, ax = plt.subplots()

    for label, file in files:
        path = os.path.join(dir_res, file)
        hvs = np.load(path)
        i, v = hvs.T
        i = np.arange(len(v)) # Generationsに変更
        ax.plot(i, v, label=label)

    ax.legend()
    # ax.set_xlabel('Evaluation Step')
    ax.set_xlabel('Number of Generations')
    ax.set_ylabel('Hypervolume')

    if show:
        plt.show()
    else:
        basename = os.path.splitext(os.path.basename(file))[0]
        plt.savefig(os.path.join(out, 'case0_hv.png'))


def plot_obj(out='result', iobj=0, show=False):
    # iobj = 0 # dv, RCS
    name = ('dv', 'rcs')[iobj]

    dir_res = 'result'
    files = [
        ('NSGA-II', 'nsga2/pop_size500_epoch100_20181121_1010_data.npy'),
        ('MOEA/D case1', 'moead/pop_size500_epoch100_20181121_2320_data.npy'), # all
        ('MOEA/D case2', 'moead/pop_size500_epoch100_20181121_1006_data.npy'), # new
        # ('MOEA/D case1', 'moead/pop_size500_epoch100_20181122_0905_data.npy'), # all(scale)
        # ('MOEA/D case2', 'moead/pop_size500_epoch100_20181122_0909_data.npy'), # new(scale)
    ]

    fig, ax = plt.subplots()

    for label, file in files:
        path = os.path.join(dir_res, file)
        data = np.load(path)
        i, *obj = data.T
        i = np.arange(len(obj[iobj])) # Generationsに変更
        ax.plot(i, obj[iobj], label=label)

    ax.legend(prop={'size': 22})
    # ax.set_xlabel('Evaluation Step')
    ax.set_xlabel('Number of Generations')
    # ax.set_ylabel('Hypervolume')
    if iobj == 0:
        ax.set_ylabel('Minimum value of $\\Delta\\it{V}$, $km/s$')
    else:
        ax.set_ylabel('Maximum value of $RCS$, $m^2$')

    # ax.set_xlabel('$\\Delta\\it{V}$, $km/s$')
    # ax.set_ylabel('$RCS$, $m^2$')

    if show:
        plt.show()
    else:
        basename = os.path.splitext(os.path.basename(file))[0]
        plt.savefig(os.path.join(out, f'case0_obj_{name}.png'))


def plot_ndom(out='result', show=False):
    dir_res = 'result'
    files = [
        ('NSGA-II', 'nsga2/pop_size500_epoch100_20181121_1010_ndom.npy'),
        ('MOEA/D case1', 'moead/pop_size500_epoch100_20181121_2320_ndom.npy'), # all
        ('MOEA/D case2', 'moead/pop_size500_epoch100_20181121_1006_ndom.npy'), # new
    ]

    fig, ax = plt.subplots()

    for label, file in files:
        path = os.path.join(dir_res, file)
        data = np.load(path)
        i, ndom = data.T
        i = np.arange(len(ndom)) # Generationsに変更
        ax.plot(i, ndom, label=label)

    ax.legend()
    ax.set_xlabel('Number of Generations')
    ax.set_ylabel('Number of Non-Dominated Solutions')

    if show:
        plt.show()
    else:
        basename = os.path.splitext(os.path.basename(file))[0]
        plt.savefig(os.path.join(out, f'case0_ndom.png'))


def __test__():
    pass


def get_args():
    '''
    docstring for get_args.
    '''
    parser = argparse.ArgumentParser()
    # parser.add_argument('out', nargs='?', default='new_script',
    #                     help='Filename of the new script')
    parser.add_argument('mode', nargs='?', default='new_script',
                        help='Filename of the new script')
    parser.add_argument('--force', '-f', action='store_true',
                        help='Force')
    parser.add_argument('--test', '-t', action='store_true',
                        help='Run as test mode')
    args = parser.parse_args()
    return args


def main():
    '''
    docstring for main.
    '''
    args = get_args()

    if args.test:
        __test__()
        return

    case.pltenv()
    if args.mode == 'hv':
        plot_hv()
    elif args.mode == 'o0':
        plot_obj(iobj=0)
    elif args.mode == 'o1':
        plot_obj(iobj=1)
    elif args.mode == 'nd':
        plot_ndom()


if __name__ == '__main__':
    main()
