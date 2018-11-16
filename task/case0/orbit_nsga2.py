import argparse
import glob
import os
import shutil
import sys
from operator import attrgetter

import numpy as np
import matplotlib.pyplot as plt

from eclib.benchmarks import rosenbrock, zdt1, zdt2, zdt3, zdt4, zdt6
from eclib.operations import UniformInitializer
from eclib.operations import RouletteSelection
from eclib.operations import TournamentSelection
from eclib.operations import TournamentSelectionStrict
from eclib.operations import TournamentSelectionDCD
from eclib.operations import BlendCrossover
from eclib.operations import SimulatedBinaryCrossover
from eclib.operations import PolynomialMutation
from eclib.optimizers import NSGA2
from eclib.base import Individual

import myutils as ut
import case0 as orbitlib


class Genome(tuple):
    ''' 進化計算遺伝子
    '''
    def __init__(self, iterable):
        self._state = 1

    def get_gene(self):
        self._state = 1 - self._state
        return self[self._state]

    def get_ivalue(self):
        return np.array([0, *self[0]])

    def get_dvalue(self):
        return self.weight * self[1]

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
        iarray = np.random.choice(range(self._imax), self._size) + 1
        darray = np.random.rand(self._size * 2)
        return Genome((iarray, darray))


class Crossover(object):
    def __init__(self, ix, dx):
        self._ix = ix
        self._dx = dx

    def __call__(self, origin):
        x1 = origin[0].get_gene() # Genome
        x2 = origin[1].get_gene() # Genome

        yi1, yi2 = self._ix((x1, x2)) # ndarray(int)
        yd1, yd2 = self._dx((x1, x2)) # ndarray(double)

        return Genome((yi1, yd1)), Genome((yi2, yd2))


class Mutation(object):
    def __init__(self, im, dm):
        self._im = im
        self._dm = dm

    def __call__(self, genome):
        yi = self._im(genome.get_gene()) # ndarray(int)
        yd = self._dm(genome.get_gene()) # ndarray(double)

        return Genome((yi, yd))


class Problem(object):
    def __init__(self, size):
        # CDLLインスタンス作成
        libname = 'case0.dll'
        loader_path = '.'
        cdll = np.ctypeslib.load_library(libname, loader_path)

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

        self._size = size
        self.function = f_call_problem
        self.imax = n_debris

    def __call__(self, genome):
        # 関数呼び出し
        order = genome.get_ivalue()
        params = genome.get_dvalue()
        delv, rcs = self.function(order, params)
        return delv, rcs

    def __reduce_ex__(self, protocol):
        return type(self), (self._size,)


################################################################################

def identity(x):
    return x

def istub(x):
    return x[0].get_gene(), x[1].get_gene()

def clip(x):
    return np.clip(x, 0.0, 1.0)


################################################################################

class NSGA_ENV(object):
    def __init__(self, epoch):
        # パラメータ
        n_dim = 5
        pop_size = 100
        epoch = 250

        # 問題
        problem = Problem(n_dim)
        initializer = Initializer(n_dim, problem.imax)

        # selection = TournamentSelection(ksize=2)
        selection = TournamentSelectionStrict(ksize=2)

        # selection = TournamentSelectionDCD()
        # crossover = BlendCrossover(alpha=0.5)

        crossover = Crossover(ix=orbitlib.OrderCrossover(),
                              dx=SimulatedBinaryCrossover(rate=0.9, eta=20))

        mutation = Mutation(im=identity,
                            dm=PolynomialMutation(rate=1/n_dim, eta=20))

        optimizer = NSGA2(pop_size, selection, crossover, mutation)
        optimizer.setup(problem)

        ### Additional setting ###
        Genome.set_weight([10, 3600] * n_dim)
        Individual.set_weight([1, -1])
        optimizer.initializer = initializer
        optimizer.n_cycle = pop_size // 2
        # optimizer.alternation = 'replace'
        ##########################

        self.optimizer = optimizer

    def __enter__(self):
        return self.optimizer

    def __exit__(self, exc_type, exc_value, traceback):
        self.optimizer.clear()


def ga_main01(out='result', clear_directory=False):
    ''' GAテスト & プロット
    '''
    if clear_directory and os.path.isdir(out):
        shutil.rmtree(out)

    epoch = 200
    save_trigger = lambda i: i == epoch # 最後だけ

    with NSGA_ENV(epoch=epoch) as optimizer:
        with ut.stopwatch('main'):
            # GA開始
            # 初期集団生成
            optimizer.create_initial_population()

            # 進化
            for i in range(1, epoch + 1):
                optimizer.advance()
                print('epoch:', i, 'popsize:', len(optimizer.population),
                      end='\r')
                if save_trigger(i):
                    optimizer.save(file=os.path.join(out, f'epoch{i}.pickle'))


        # elite = optimizer.get_elite()
        # history = optimizer.history
        # def best(pop):
        #     return [x for x in pop if x.rank == 1][0]()
        # bests = np.array([best(pop) for pop in history])

        # first_population = optimizer[0]

        last_population = optimizer.get_individuals()
        last_population.sort(key=lambda x: x.value)
        # optimal_front = get_optomal_front('pareto_front/zdt1_front.json')

        ### TEMP: check stat ###
        # print("Convergence: ", convergence(last_population, optimal_front))
        # print("Diversity: ", diversity(last_population, optimal_front[0], optimal_front[-1]))
        ########################

        ### TEMP: plot front ###
        x, y = np.array([x.value for x in last_population]).T

        # plt.scatter(optimal_front[:, 0], optimal_front[:, 1], c='r')

        plt.scatter(x, y, c='b')
        plt.axis("tight")
        # plt.xlim((0, 1))
        # plt.ylim((0, 1))
        plt.show()
        ########################


def ga_main02(out='result', clear_directory=False):
    ''' GAテスト & プロット
    '''
    if clear_directory and os.path.isdir(out):
        shutil.rmtree(out)

    epoch = 250
    save_trigger = lambda i: i == epoch # 最後だけ
    optimal_front = get_optomal_front()
    stat = []

    with NSGA_ENV(epoch) as optimizer:
        for rep in range(100):
            with ut.stopwatch(f'epoch{epoch+1}'):
                optimizer.create_initial_population()
                for i in range(1, epoch + 1):
                    optimizer.advance()
                    print('epoch:', i, 'popsize:', len(optimizer.population), end='\r')

            last_population = optimizer.get_individuals()
            last_population.sort(key=lambda x: x.value)

            conv = convergence(last_population, optimal_front)
            div = diversity(last_population, optimal_front[0], optimal_front[-1])
            stat.append((conv, div))

            print("Convergence: ", conv)
            print("Diversity: ", div)

    print('=' * 20, 'Average', '=' * 20)
    print("Convergence: ", np.mean([x[0] for x in stat]))
    print("Diversity: ",  np.mean([x[1] for x in stat]))


def ga_main1(out='result'):
    file = ut.fsort(glob.glob(os.path.join(out, f'epoch*.pickle')))[-1]
    optimizer = NSGA2.load(file=file)
    print('epoch:', len(optimizer))

    elite = optimizer.get_elite()
    print('elite:', len(elite))

    for epoch in range(0, len(optimizer)-1, 10):
        print(epoch, end='\n')
        plt.cla()
        population = optimizer[epoch]
        population = optimizer.calc_rank(population)

        for i in range(1, 2):
            front = [x for x in population if x.rank == i]
            # print(sorted([x.rank for x in population]))
            if not front:
                continue
            x, y = np.array([x.data.value for x in front]).T
            # print(front)
            # return
            plt.scatter(x, y, label=f'{front[0][0]:.3f}')

        plt.legend()
        plt.xlim((0, 50))
        if epoch < len(optimizer)-1:
            plt.pause(0.2)
        else:
            plt.show()


def ga_main2(out='result'):
    file = ut.fsort(glob.glob(os.path.join(out, f'epoch*.pickle')))[-1]
    optimizer = NSGA2.load(file=file)

    population = optimizer[-1]
    front = [x for x in population if x.rank == 1]
    front.sort(key=attrgetter('data.value'))

    for ind in front:
        print(ind.value, ind.data.value)

    crowdings = [ind.value[1] for ind in front]

    fig, axes = plt.subplots(2)
    axes[0].plot(crowdings)

    x, y = np.array([x.data.value for x in front]).T
    im = axes[1].scatter(x, y, c=crowdings, cmap='jet')
    # plt.xlim((0, 1))
    # plt.ylim((0, 1))
    fig.colorbar(im)
    plt.show()

    # for ind in first_population:
    #     print(ind(), ind.rank, ind.fitness)
    # return

    # last_population = optimizer.population
    # x, y = np.array([x() for x in last_population]).T
    # plt.scatter(x, y)

    # plt.xlim((0, 1))
    # plt.ylim((0, 1))
    # plt.show()


################################################################################
### TEMP ###
################################################################################

def get_optomal_front(file):
    import json
    with open(file, 'r') as optimal_front_data:
        optimal_front = json.load(optimal_front_data)
    # Use 500 of the 1000 points in the json file
    optimal_front = sorted(optimal_front[i]
                           for i in range(0, len(optimal_front), 2))
    return np.array(optimal_front)


def diversity(first_front, first, last):
    from math import hypot
    df = hypot(first_front[0][0] - first[0],
               first_front[0][1] - first[1])
    dl = hypot(first_front[-1][0] - last[0],
               first_front[-1][1] - last[1])
    dt = [hypot(f[0] - s[0],
                f[1] - s[1])
          for f, s in zip(first_front[:-1], first_front[1:])]

    if len(first_front) == 1:
        return df + dl

    dm = sum(dt)/len(dt)
    di = sum(abs(d_i - dm) for d_i in dt)
    delta = (df + dl + di)/(df + dl + len(dt) * dm )
    return delta


def convergence(first_front, optimal_front):
    from math import sqrt
    distances = []

    for ind in first_front:
        distances.append(float("inf"))
        for opt_ind in optimal_front:
            dist = 0.
            for i in range(len(opt_ind)):
                dist += (ind[i] - opt_ind[i])**2
            if dist < distances[-1]:
                distances[-1] = dist
        distances[-1] = sqrt(distances[-1])

    return sum(distances) / len(distances)


################################################################################

def __test__():
    genome = Genome((0, 1))
    print(genome.get_gene())
    print(genome.get_gene())


def get_args():
    '''
    docstring for get_args.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('method', nargs='?', default='',
                        help='Main method type')
    parser.add_argument('--out', '-o', default='result',
                        help='Filename of the new script')
    parser.add_argument('--clear', '-c', action='store_true',
                        help='Remove output directory before start')
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

    # print(sys.getrecursionlimit())
    sys.setrecursionlimit(10000)

    args = get_args()
    out = args.out
    clear = args.clear

    if args.test:
        __test__()
        return

    if args.method == '0':
        ga_main(out=out, clear_directory=clear)
    elif args.method == '01':
        ga_main01(out=out, clear_directory=clear)
    elif args.method == '02':
        ga_main02(out=out, clear_directory=clear)
    elif args.method == '1':
        ga_main1()
    elif args.method == '2':
        ga_main2()


if __name__ == '__main__':
    main()
