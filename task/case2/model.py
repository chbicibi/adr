import argparse
import ctypes
import glob
import os
import shutil
import sys
from itertools import chain
from operator import attrgetter, itemgetter

import numpy as np

from eclib.base import Individual
from eclib.base import Environment
from eclib.base import Creator
from eclib.operations import UniformInitializer
from eclib.operations import TournamentSelection
from eclib.operations import SimulatedBinaryCrossover
from eclib.operations import PolynomialMutation
from eclib.optimizers import NSGA2
from eclib.optimizers import MOEAD

# テスト関数
from eclib.benchmarks import rosenbrock, zdt1, zdt2, zdt3, zdt4, zdt6

import utils as ut
import problem as orbitlib


################################################################################

class Problem(object):
    def __init__(self, size=0):
        f_call_problem, n_debris = get_fortran_function()
        self._size = size
        self.function = f_call_problem
        self.imax = n_debris

    def __call__(self, params):
        # 関数呼び出し
        indexes = [1, 10]
        delv, delt = self.function(indexes, params)
        return delv, delt

    def __reduce_ex__(self, protocol):
        return type(self), (self._size,)


def get_fortran_function(cache=[]):
    if not cache:
        # CDLLインスタンス作成
        libname = 'problem.dll'
        loader_path = 'fortran'
        cdll = np.ctypeslib.load_library(libname, loader_path)
        # cdll = ctypes.WinDLL(libname)

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
        cache.extend([f_call_problem, n_debris])
    else:
        f_call_problem, n_debris = cache

    return f_call_problem, n_debris


def manual_problem(params):
    function, _ = get_fortran_function()
    indexes = [1, 2]
    delv, delt = function(indexes, params)
    return delv, delt



################################################################################

class Optimize_ENV(object):
    def __init__(self, method, popsize=100, **kwargs):
        method = method.lower()
        if 'nsga' in method:
            opt_cls = NSGA2
            ga_ops = {}
        elif 'moea' in method:
            opt_cls = MOEAD
            if kwargs['ksize']:
                ksize = kwargs['ksize']
            else:
                ksize = 5
            ga_ops = {'ksize': ksize}
        else:
            print('Unknown method name:', method)
            raise RuntimeError

        ########################################################################
        # 最適化問題設定
        ########################################################################
        # 設計範囲
        # [出発日, 軌道長半径比, ランデブー中心角, 軌道面偏向角比1, 軌道面偏向角比2, 周回数1, 周回数2]
        low_bounds = [0, 0.5, 0.3, 0.1, 0.1, 0.0, 0.0]
        upp_bounds = [30, 1.5, 0.7, 0.8, 0.8, 10.0, 10.0]

        # 設計変数の次元数
        n_dim = len(low_bounds)

        # 最適化重み(正=>最小化, 負=>最大化)
        opt_weight = [1, 1]

        # 問題関数
        # problem = Problem()
        problem = manual_problem # 上の行と同じ

        ########################################################################
        # 最適化アルゴリズム構築
        ########################################################################
        with Environment() as env:
            # 個体クラス
            indiv_pool = env.register(Individual)

            # 遺伝子初期化クラス
            initializer = UniformInitializer(n_dim)

            # GAオペレータ指定
            ga_ops = {
              'selection': TournamentSelection(ksize=2),
              'crossover': SimulatedBinaryCrossover(rate=0.9, eta=20),
              'mutation': PolynomialMutation(rate=0.1, eta=20),
              **ga_ops
            }

            # GAクラス
            optimizer = opt_cls(problem=problem, pool=indiv_pool, **ga_ops)

            ### Additional setting ###
            # 設計範囲設定
            indiv_pool.cls.set_bounds(low_bounds, upp_bounds) # (下限，上限)

            # 最適化重み付け(正=>最小化, 負=>最大化)
            indiv_pool.cls.set_weight(opt_weight)

            # 親個体選択時の初期化周期
            optimizer.n_cycle = popsize // 2
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


################################################################################

def __test__():
    pass


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
    args = get_args()

    if args.test:
        __test__()
        return


if __name__ == '__main__':
    main()
