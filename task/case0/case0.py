#! /usr/bin/env python3

'''
1. 初期集団生成
2.
'''

import argparse
import os
import shutil

import numpy as np
import matplotlib.pyplot as plt

from benchmarks import zdt1, rosenbrock
from nsga2 import NSGA2
from operations.selection import roulette
from operations.crossover import blx

# import operations.crossover


class Initializer(object):
    ''' [0, 1)の範囲の一様乱数による実数配列を返す '''

    def __init__(self, size):
        self._size = size

    def __call__(self):
        return np.random.rand(self._size)


# class Individual(object):
#     ''' 進化計算個体
#     saku
#     '''

#     def __init__(self):
#         self.genotype

#     def __call__(self):
#         return np.random.rand(self._size)


# class Population(object):
#     ''' 進化計算集団 '''

#     def __init__(self, n):
#         self._size = n

#     def __call__(self):
#         return np.random.rand(self._size)


################################################################################

def identity(x):
    return x


def clip(x):
    return np.clip(x, 0, 1)


################################################################################

def ga_main():
    # パラメータ
    pop_size = 100
    epoch = 50

    # 問題
    problem = zdt1
    selection = roulette
    crossover = blx
    mutation = clip

    optimizer = NSGA2(pop_size, selection, crossover, mutation)
    optimizer.set_initializer(Initializer(3))
    optimizer.setup(problem)

    # GA開始
    # 初期集団生成
    optimizer.start()

    # 進化
    for i in range(1, epoch + 1):
        optimizer.advance()
        if i == 1 or i % 1 == 0:
            print('epoch:', i)
            optimizer.save(file=f'result/epoch{i}.pickle')

    # elite = optimizer.get_elite()
    # history = optimizer.history
    # def best(pop):
    #     return [x for x in pop if x.rank == 1][0]()
    # bests = np.array([best(pop) for pop in history])

    # first_population = optimizer[0]

    last_population = optimizer.current_population
    x, y = np.array([x() for x in last_population]).T
    plt.scatter(x, y)
    # plt.xlim((0, 1))
    # plt.ylim((0, 1))
    plt.show()


def ga_main1():
    optimizer = NSGA2.load(file=f'result/epoch{50}.pickle')

    elite = optimizer.get_elite()
    print('elite:', len(elite))

    for epoch in range(50):
        plt.cla()
        first_population = optimizer[epoch]

        for i in range(10):
            pop = [x for x in first_population if x.rank == i]
            if not pop:
                continue
            x, y = np.array([x() for x in pop]).T
            plt.scatter(x, y, label=f'{pop[0].fitness:.3f}')

        plt.legend()
        plt.pause(1)

    # for ind in first_population:
    #     print(ind(), ind.rank, ind.fitness)
    # return

    # last_population = optimizer.current_population
    # x, y = np.array([x() for x in last_population]).T
    # plt.scatter(x, y)

    # plt.xlim((0, 1))
    # plt.ylim((0, 1))
    # plt.legend()
    # plt.show()


################################################################################

def __test__():
    from scipy.interpolate import RegularGridInterpolator
    from mpl_toolkits.mplot3d import Axes3D

    x = np.arange(0, 1, 0.2)
    y = np.arange(0, 1, 0.2)

    X, Y = np.meshgrid(x, y)
    Z = np.array([[rosenbrock([vx, vy]) for vy in y] for vx in x])

    # print(X.shape)
    # # print(y.shape)
    # print(Z.shape)
    # return
    xx = np.arange(0, 0.8, 0.01)
    yy = np.arange(0, 0.8, 0.01)
    XX, YY = np.meshgrid(xx, yy)
    # ax.scatter3D(*map(lambda a: a.reshape(-1), (X, Y, Z)))

    zz = RegularGridInterpolator((x, y), Z)

    ZZ = np.array([[zz([vx, vy]) for vy in yy] for vx in xx])[:, :, 0]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # ax.plot_surface(XX, YY, ZZ)
    ax.scatter3D(*map(lambda a: a.reshape(-1), (XX, YY, ZZ)))

    ax.set_title("Scatter Plot")
    plt.show()

    # plt.scatter(z[:, 0], z[:, 1])
    # plt.xlim((0, 1))
    # # plt.ylim((0, 1))
    # plt.show()


# def __test__():
#     a = list(range(5))
#     print(a.pop(3))
#     print(a)


def get_args():
    '''
    docstring for get_args.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('method', nargs='?', default='',
                        help='Main method type')
    parser.add_argument('--out', '-o', default='new_script',
                        help='Filename of the new script')
    parser.add_argument('--force', '-f', action='store_true',
                        help='force')
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

    if args.method == '0':
        ga_main()
    elif args.method == '1':
        ga_main1()


if __name__ == '__main__':
    main()
