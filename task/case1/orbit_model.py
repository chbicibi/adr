import argparse
import ctypes
import glob
import os
import shutil
import sys
from itertools import chain
from operator import attrgetter, itemgetter

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as plc
from matplotlib import animation as anim
from matplotlib import rcParams

from eclib.benchmarks import rosenbrock, zdt1, zdt2, zdt3, zdt4, zdt6
from eclib.operations import UniformInitializer
from eclib.operations import RandomSelection
from eclib.operations import RouletteSelection
from eclib.operations import TournamentSelection
from eclib.operations import TournamentSelectionStrict
from eclib.operations import TournamentSelectionDCD
from eclib.operations import BlendCrossover
from eclib.operations import SimulatedBinaryCrossover
from eclib.operations import OrderCrossover
from eclib.operations import PolynomialMutation
from eclib.operations import SwapMutation
from eclib.optimizers import NSGA2
from eclib.optimizers import MOEAD
from eclib.base import Individual
from eclib.base import NoisingIndividual
from eclib.base import Environment
from eclib.base import Creator
from eclib.base import Population
from eclib.base import HyperVolume

import utils as ut
import problem as orbitlib


def init_ffmpeg():
    rcParams["animation.ffmpeg_path"] = "ffmpeg.exe"
    FFMpegWriter = anim.writers['ffmpeg']
    return FFMpegWriter


def show_anim(fig, update, frames=1000, init_func=lambda:None, interval=8, file=None, fps=2):
  ani = anim.FuncAnimation(fig, update, frames=frames, init_func=init_func, interval=interval)
  if file:
    ani.save(file, writer=FFMpegWriter(fps=fps))
  else:
    plt.show()


################################################################################

class Genome(tuple):
    ''' 進化計算遺伝子
    '''
    def __init__(self, iterable):
        self._state = 1

    def get_gene(self, dtype=None):
        if dtype == 'i':
            return self[0]
        elif dtype == 'd':
            return self[1]
        self._state = 1 - self._state
        return self[self._state]

    def get_ivalue(self):
        return np.array([0, *self[0]])

    def get_dvalue(self):
        # print('get_dvalue')
        # print(self.weight)
        # print(self[1])
        # print('#')
        return self.weight * self[1]

    def set_state(self, state):
        self._state = state

    @classmethod
    def set_weight(cls, weight):
        cls.weight = np.array(weight)


class Initializer(object):
    ''' [0, 1)の範囲の一様乱数による実数配列を返す
    '''
    def __init__(self, size, imax):
        self._size = size
        self._imax = imax

    def __call__(self):
        iarray = np.random.choice(range(self._imax), self._size,
                                  replace=False) + 1
        darray = np.random.rand(self._size * 2)
        return Genome((iarray, darray))


class Crossover(object):
    def __init__(self, ix, dx):
        self._ix = ix
        self._dx = dx

    def __call__(self, origin):
        ''' origin.__class__ => [Genome] '''
        xi1 = origin[0].get_gene('i') # Genome
        xd1 = origin[0].get_gene('d') # Genome
        xi2 = origin[1].get_gene('i') # Genome
        xd2 = origin[1].get_gene('d') # Genome
        # exit()
        # print(xi1, xd1)
        # exit()

        # yi1, yi2 = self._ix((xi1, xi2)) # ndarray(int)
        yi1, yi2 = xi1, xi2 # ndarray(int)
        yd1, yd2 = self._dx((xd1, xd2)) # ndarray(double)

        return Genome((yi1, yd1)), Genome((yi2, yd2))


class Mutation(object):
    def __init__(self, im, dm):
        self._im = im
        self._dm = dm

    def __call__(self, genome):
        yi = self._im(genome.get_gene('i')) # ndarray(int)
        yd = self._dm(genome.get_gene('d')) # ndarray(double)

        return Genome((yi, yd))

FUNC_FORTRAN = None
class Problem(object):
    def __init__(self, size, arg):
        global FUNC_FORTRAN

        if not FUNC_FORTRAN:
            # CDLLインスタンス作成
            libname = 'case0.dll'
            loader_path = '.'
            # cdll = np.ctypeslib.load_library(libname, loader_path)
            cdll = ctypes.WinDLL(libname)

            # 関数取得
            f_initialize = orbitlib.get_f_initialize(cdll)
            f_init_debri = orbitlib.get_f_init_debri(cdll)
            f_call_problem = orbitlib.get_f_call_problem(cdll)

            # 初期化: 開始時刻設定
            f_initialize()

            # 初期化: デブリデータ読み込み
            tle_file = '../data/debri_elements.txt'
            rcs_file = '../data/RCS_list.txt'
            n_debris = f_init_debri(tle_file, rcs_file)
            print('n_debris:', n_debris)
            FUNC_FORTRAN = (f_call_problem, n_debris)
        else:
            f_call_problem, n_debris = FUNC_FORTRAN

        self._size = size
        self.function = f_call_problem
        self.imax = n_debris
        self._arg = arg

    def __call__(self, genome):
        # 関数呼び出し
        order = genome.get_ivalue()
        params = genome.get_dvalue()
        delv, rcs = self.function(order, params)
        return delv, rcs

    def __reduce_ex__(self, protocol):
        return type(self), (self._size, self._arg)


################################################################################

def identity(x):
    return x

def istub(x):
    return x[0].get_gene(), x[1].get_gene()

def clip(x):
    return np.clip(x, 0.0, 1.0)


################################################################################

class NSGA2_ENV(object):
    def __init__(self, popsize=100, **kwargs):
        # パラメータ
        n_dim = 10
        # popsize = 200
        # epoch = 100
        imax = 100 # 対象デブリ数

        # 問題
        problem = Problem(n_dim)

        with Environment() as env:
            # 個体クラス
            indiv_type = Individual
            # 初期個体生成クラス
            indiv_pool = env.register(indiv_type)

            # 遺伝子生成クラス
            initializer = Initializer(n_dim, imax)

            # GAオペレータ指定
            # selection = TournamentSelection(ksize=2)
            selection = TournamentSelectionStrict(ksize=2)

            # selection = TournamentSelectionDCD()
            # crossover = BlendCrossover(alpha=0.5)

            crossover = Crossover(ix=OrderCrossover(rate=0.9),
                                  dx=SimulatedBinaryCrossover(rate=0.9, eta=20))

            mutation = Mutation(im=SwapMutation(rate=0.1),
                                dm=PolynomialMutation(rate=0.1, eta=20))

            ga_ops = {'selection': selection,
                      'crossover': crossover,
                      'mutation': mutation}

            optimizer = NSGA2(problem=problem, pool=indiv_pool, **ga_ops)
            # optimizer = MOEAD(problem=problem, pool=indiv_pool, ksize=5, **ga_ops)

            # optimizer = NSGA2(popsize, selection, crossover, mutation,
            #                   indiv_type=Individual)
            # optimizer.setup(problem)

            ### Additional setting ###
            Genome.set_weight([1, 3600] * n_dim)
            indiv_type.set_weight([1, -1]) # [dv, rcs]
            # optimizer.initializer = initializer
            optimizer.n_cycle = popsize // 2
            # optimizer.alternation = 'replace'
            ##########################

            # 個体生成器
            creator = Creator(initializer, indiv_pool)

            # 登録
            env.optimizer = optimizer
            env.creator = creator
            self.env = env

    def __enter__(self):
        return self.env

    def __exit__(self, exc_type, exc_value, traceback):
        pass
        # self.optimizer.clear()


class MOEAD_ENV(object):
    def __init__(self, popsize=100, ksize=5):
        # パラメータ
        n_dim = 10
        # popsize = 200
        # epoch = 100
        imax = 100 # 対象デブリ数
        # ksize = 5

        # 問題
        problem = Problem(n_dim)
        with Environment() as env:
            # 個体クラス
            indiv_type = Individual
            # 初期個体生成クラス
            indiv_pool = env.register(indiv_type)

            # 遺伝子生成クラス
            initializer = Initializer(n_dim, imax)

            # GAオペレータ指定
            # selection = TournamentSelection(ksize=2)
            selection = TournamentSelectionStrict(ksize=2)

            # selection = TournamentSelectionDCD()
            # crossover = BlendCrossover(alpha=0.5)

            crossover = Crossover(ix=OrderCrossover(rate=0.9),
                                  dx=SimulatedBinaryCrossover(rate=0.9, eta=20))

            mutation = Mutation(im=SwapMutation(rate=0.1),
                                dm=PolynomialMutation(rate=0.1, eta=20))

            ga_ops = {'selection': selection,
                      'crossover': crossover,
                      'mutation': mutation}

            # optimizer = NSGA2(problem=problem, pool=indiv_pool, **ga_ops)
            optimizer = MOEAD(problem=problem, pool=indiv_pool, ksize=5, **ga_ops)


            ### Additional setting ###
            Genome.set_weight([1, 3600] * n_dim)
            Individual.set_weight([1, -1])
            # optimizer.initializer = initializer
            optimizer.n_cycle = popsize // 2
            optimizer.alternation = 'all'
            ##########################

            # 個体生成器
            creator = Creator(initializer, indiv_pool)

            # 登録
            env.optimizer = optimizer
            env.creator = creator
            self.env = env

        # self.optimizer = optimizer

    def __enter__(self):
        return self.env

    def __exit__(self, exc_type, exc_value, traceback):
        pass
        # self.optimizer.clear()


################################################################################

def __test__():
    libname = 'case0.dll'
    loader_path = '.'
    # cdll = np.ctypeslib.load_library(libname, loader_path)
    cdll = ctypes.WinDLL(libname)

    f_initialize = orbitlib.get_f_initialize(cdll)
    f_init_debri = orbitlib.get_f_init_debri(cdll)
    f_call_problem = orbitlib.get_f_call_problem(cdll)
    f_initialize()
    print(cdll)


def get_args():
    '''
    docstring for get_args.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', '-t', action='store_true',
                        help='Run as test mode')
    args = parser.parse_args()
    return args


def main():

    if args.test:
        __test__()
        return


if __name__ == '__main__':
    main()
