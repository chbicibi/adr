import argparse
import glob
import os
import shutil
import sys
from operator import attrgetter, itemgetter

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as plc
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
from eclib.base import HyperVolume

import myutils as ut
import case0 as orbitlib


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
        x1 = origin[0].get_gene() # Genome
        x2 = origin[1].get_gene() # Genome
        # exit()

        yi1, yi2 = self._ix((x1, x2)) # ndarray(int)
        yd1, yd2 = self._dx((x1, x2)) # ndarray(double)

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
    def __init__(self, size):
        global FUNC_FORTRAN

        if not FUNC_FORTRAN:
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
            FUNC_FORTRAN = (f_call_problem, n_debris)
        else:
            f_call_problem, n_debris = FUNC_FORTRAN

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

class NSGA2_ENV(object):
    def __init__(self, pop_size=100, **kwargs):
        # パラメータ
        n_dim = 10
        # pop_size = 200
        # epoch = 100
        imax = 100 # 対象デブリ数

        # 問題
        problem = Problem(n_dim)
        initializer = Initializer(n_dim, imax)

        # selection = TournamentSelection(ksize=2)
        selection = TournamentSelectionStrict(ksize=2)

        # selection = TournamentSelectionDCD()
        # crossover = BlendCrossover(alpha=0.5)

        crossover = Crossover(ix=OrderCrossover(rate=0.9),
                              dx=SimulatedBinaryCrossover(rate=0.9, eta=20))

        mutation = Mutation(im=SwapMutation(rate=0.1),
                            dm=PolynomialMutation(rate=0.1, eta=20))

        optimizer = NSGA2(pop_size, selection, crossover, mutation,
                          indiv_type=Individual)
        optimizer.setup(problem)

        ### Additional setting ###
        Genome.set_weight([1, 3600] * n_dim)
        Individual.set_weight([1, -1]) # [dv, rcs]
        optimizer.initializer = initializer
        optimizer.n_cycle = pop_size // 2
        # optimizer.alternation = 'replace'
        ##########################

        self.optimizer = optimizer

    def __enter__(self):
        return self.optimizer

    def __exit__(self, exc_type, exc_value, traceback):
        self.optimizer.clear()


class MOEAD_ENV(object):
    def __init__(self, pop_size=100, ksize=5):
        # パラメータ
        n_dim = 10
        # pop_size = 200
        # epoch = 100
        imax = 100 # 対象デブリ数
        # ksize = 5

        # 問題
        problem = Problem(n_dim)
        initializer = Initializer(n_dim, imax)

        # selection = TournamentSelection(ksize=2)
        # selection = TournamentSelectionStrict(ksize=2)
        selection = RandomSelection()

        # selection = TournamentSelectionDCD()
        # crossover = BlendCrossover(alpha=0.5)

        crossover = Crossover(ix=OrderCrossover(rate=0.9),
                              dx=SimulatedBinaryCrossover(rate=0.9, eta=20))

        mutation = Mutation(im=SwapMutation(rate=0.1),
                            dm=PolynomialMutation(rate=0.1, eta=20))

        optimizer = MOEAD(pop_size, selection, crossover, mutation, ksize)
        optimizer.setup(problem)

        ### Additional setting ###
        Genome.set_weight([1, 7200] * n_dim)
        Individual.set_weight([1, -1])
        optimizer.initializer = initializer
        optimizer.n_cycle = pop_size // 2
        optimizer.alternation = 'new'
        ##########################

        self.optimizer = optimizer

    def __enter__(self):
        return self.optimizer

    def __exit__(self, exc_type, exc_value, traceback):
        self.optimizer.clear()


################################################################################

def ga_main1(model, out='result', clear_directory=False):
    ''' GAテスト & プロット
    '''
    if clear_directory and os.path.isdir(out):
        shutil.rmtree(out)


    pop_size = 500
    epoch = 100
    ksize = 5
    save_trigger = lambda i: i == epoch # 最後だけ

    model_env = {'nsga2':NSGA2_ENV, 'moead':MOEAD_ENV}[model]

    with model_env(pop_size=pop_size, ksize=ksize) as optimizer:
        with ut.stopwatch('main'):
            # GA開始
            # 初期集団生成
            optimizer.init_population()

            # 進化
            for i in range(1, epoch + 1):
                optimizer.advance()
                print('epoch:', i, 'popsize:', len(optimizer.population),
                      end='\r')
                if save_trigger(i):
                    file = f'pop_size{pop_size}_epoch{i}_{ut.strnow("%Y%m%d_%H%M")}.pickle'
                    optimizer.save(file=os.path.join(out, file))

    ga_result11(model, out=out, show=False)


def ga_main2(out='result', clear_directory=False):
    ''' GAテスト & プロット
    '''
    if clear_directory and os.path.isdir(out):
        shutil.rmtree(out)

    epoch = 250
    save_trigger = lambda i: i == epoch # 最後だけ
    optimal_front = get_optomal_front()
    stat = []

    with MOEAD_ENV() as optimizer:
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


def ga_result1(out='result'):
    ''' アニメーション表示
    '''
    file = ut.fsort(glob.glob(os.path.join(out, f'epoch*.pickle')))[-1]
    optimizer = MOEAD.load(file=file)
    print('epoch:', len(optimizer))

    elite = optimizer.get_elite()
    print('elite:', len(elite))

    imax = len(optimizer) - 1

    for epoch in range(0, imax, 20):
        print(epoch, end='\n')
        # plt.cla()
        population = optimizer[epoch]
        # population = optimizer.calc_rank(population)

        front = [x for x in population]
        x, y = np.array([x.data.value for x in front]).T
        plt.scatter(x, y, color='blue', alpha=epoch/imax, label=f'epoch{epoch}')

        # plt.legend()
        plt.xlim((0, 50))
        # plt.ylim((0, 50))
        plt.pause(0.2)
    plt.show()

def ga_result11(model, out='result', show=False):
    ''' 一括表示
    '''
    # モデル読み込み
    model_cls = {'nsga2':NSGA2, 'moead':MOEAD}[model]
    files = ut.fsort(glob.glob(os.path.join(out, f'*epoch*.pickle')))
    for i, file in enumerate(files):
        print(f'[{i}]', file)
    print('select file')
    n = int(input())
    if n < 0:
        return
    file = files[n]
    print('file:', file)
    optimizer = model_cls.load(file=file)

    # モデルパラメータ表示
    print('[PARAMS]')
    print('popsize:', optimizer.popsize)
    print('epoch:', len(optimizer))
    print('alternation:', optimizer.alternation)
    if hasattr(optimizer, 'ksize'):
        print('ksize:', optimizer.ksize)
    print('Indiv:', optimizer.indiv_type.__name__)
    optimizer.calc_rank(optimizer.population)
    elite = optimizer.get_elite()
    print('elite:', len(elite))

    imax = len(optimizer)
    # imax = 51
    data = []

    for epoch in range(0, imax, 1):
        print(epoch, end='\r')
        # plt.cla()
        population = optimizer[epoch]
        front = [x for x in population]
        data.extend([[*x.data.value, epoch] for x in front])

    x, y, c = np.array(data).T
    fig, ax = plt.subplots()
    # color_list = [(0, '#AADDFF'), (0.33, '#002288'), (0.66, '#EEDD22'), (1, '#CC0000')]
    color_list = [(0, '#0022AA'), (0.3, '#11AA00'), (0.5, '#EEDD22'), (1, '#CC0000')]
    color_list = [(0, 'blue'), (0.3, 'green'), (0.5, 'yellow'), (1, 'red')]
    cmap = plc.LinearSegmentedColormap.from_list('custom_cmap', color_list, N=25)
    im = ax.scatter(x, y, c=c, cmap=cmap, marker='o', alpha=1)

    print('min dV=', x.min())
    print('max RCS=', y.max())

    # plt.legend()
    ax.set_xlabel('$\\Delta\\it{V}$, $km/s$')
    ax.set_ylabel('$RCS$, $m^2$')

    ax.set_xlim((0, 100))
    ax.set_ylim((0, 0.45))

        # plt.pause(0.2)
    cbar = fig.colorbar(im, label='Generation')
    cbar.set_label('Generation', size=24)
    plt.tight_layout()

    # show or save
    if show:
        plt.show()
    else:
        basename = os.path.splitext(os.path.basename(file))[0]
        plt.savefig(os.path.join('result', f'front_{model}_{basename}.png'))


def ga_result12(model, out='result', show=False):
    ''' 最適値プロット
    '''
    # モデル読み込み
    model_cls = {'nsga2':NSGA2, 'moead':MOEAD}[model]
    files = ut.fsort(glob.glob(os.path.join(out, f'*epoch*.pickle')))
    for i, file in enumerate(files):
        print(f'[{i}]', file)
    print('select file')
    n = int(input())
    if n < 0:
        return
    file = files[n]
    print('file:', file)
    optimizer = model_cls.load(file=file)

    # モデルパラメータ表示
    print('[PARAMS]')
    print('popsize:', optimizer.popsize)
    print('epoch:', len(optimizer))
    print('alternation:', optimizer.alternation)
    print('Indiv:', optimizer.indiv_type.__name__)
    optimizer.calc_rank(optimizer.population)
    elite = optimizer.get_elite()
    print('elite:', len(elite))

    imax = len(optimizer)
    # imax = 100
    idx = 0

    cache = file.replace('.pickle', '_data.npy')
    if os.path.isfile(cache) and not force:
        data = np.load(cache)

    else:
        data = []

        for epoch in range(0, imax, 1):
            print(epoch, end='\r')
            # plt.cla()
            population = optimizer[epoch]
            indivs = [fit.get_indiv() for fit in population]
            dv = min(map(itemgetter(0), indivs))
            rcs = max(map(itemgetter(1), indivs))
            if hasattr(population[0], 'id'):
                idx = max(*map(attrgetter('id'), population), idx)
            else:
                idx = max(*map(attrgetter('id'), indivs), idx)
            data.append([idx, dv, rcs])
        np.save(cache, data)

    i, x, y = np.array(data).T
    fig, ax = plt.subplots()
    # color_list = [(0, '#AADDFF'), (0.33, '#002288'), (0.66, '#EEDD22'), (1, '#CC0000')]
    # color_list = [(0, '#0022AA'), (0.3, '#11AA00'), (0.5, '#EEDD22'), (1, '#CC0000')]
    # color_list = [(0, 'blue'), (0.3, 'green'), (0.5, 'yellow'), (1, 'red')]
    # cmap = plc.LinearSegmentedColormap.from_list('custom_cmap', color_list, N=25)
    im = ax.plot(i, x)
    im = ax.plot(i, y*100)

    print('min dV=', x.min())
    print('max RCS=', y.max())

    # plt.legend()
    # ax.set_xlabel('$\\Delta\\it{V}$, $km/s$')
    # ax.set_ylabel('$RCS$, $m^2$')

    # ax.set_xlim((0, 100))
    # ax.set_ylim((0, 0.45))

    # cbar = fig.colorbar(im, label='generation')

    if show:
        plt.show()
    else:
        basename = os.path.splitext(os.path.basename(file))[0]
        plt.savefig(os.path.join('result', f'front_{model}_{basename}.png'))


def ga_result13(model, out='result', show=False):
    ''' HVプロット
    '''
    # モデル読み込み
    model_cls = {'nsga2':NSGA2, 'moead':MOEAD}[model]
    files = ut.fsort(glob.glob(os.path.join(out, f'*epoch*.pickle')))
    for i, file in enumerate(files):
        print(f'[{i}]', file)
    print('select file')
    n = int(input())
    if n < 0:
        return
    file = files[n]
    print('file:', file)
    optimizer = model_cls.load(file=file)

    ref_point = [100, 0]
    hypervolume = HyperVolume(ref=ref_point)

    # モデルパラメータ表示
    print('[PARAMS]')
    print('popsize:', optimizer.popsize)
    print('epoch:', len(optimizer))
    print('alternation:', optimizer.alternation)
    print('Indiv:', optimizer.indiv_type.__name__)
    optimizer.calc_rank(optimizer.population)
    elite = optimizer.get_elite()
    print('elite:', len(elite))

    imax = len(optimizer)
    # imax = 100
    idx = 0
    force = False

    cache = file.replace('.pickle', '_hv.npy')
    if os.path.isfile(cache) and not force:
        hvs = np.load(cache)

    else:
        hvs = []
        population = []
        for epoch in range(0, imax, 1):
            # print(epoch, end='\r')
            # plt.cla()
            population = optimizer[epoch] + population
            optimizer.calc_rank(population)
            population = [fit for fit in population if fit.rank == 1]
            hv = hypervolume(population)
            print('epoch:', epoch, 'HV=', hv)
            if hasattr(population[0], 'id'):
                idx = max(*map(attrgetter('id'), population), idx)
            else:
                idx = max(*map(attrgetter('id'), indivs), idx)
            hvs.append([idx, hv])
        np.save(cache, hvs)

    i, v = np.array(hvs).T
    fig, ax = plt.subplots()
    im = ax.plot(i, v, label=model)

    ax.legend()
    ax.set_xlabel('Evaluation Step')
    ax.set_ylabel('Hypervolume')

    if show:
        plt.show()
    else:
        basename = os.path.splitext(os.path.basename(file))[0]
        plt.savefig(os.path.join('result', f'front_{model}_{basename}_hv.png'))


def ga_result14(model, out='result', show=False):
    ''' エリート個体数を数える
    '''
    # モデル読み込み
    model_cls = {'nsga2':NSGA2, 'moead':MOEAD}[model]
    files = ut.fsort(glob.glob(os.path.join(out, f'*epoch*.pickle')))
    for i, file in enumerate(files):
        print(f'[{i}]', file)
    print('select file')
    n = int(input())
    if n < 0:
        return
    file = files[n]
    print('file:', file)
    optimizer = model_cls.load(file=file)

    imax = len(optimizer)
    # imax = 100
    idx = 0

    cache = file.replace('.pickle', '_ndom.npy')
    if os.path.isfile(cache) and not force:
        data = np.load(cache)

    else:
        data = []

        for epoch in range(0, imax, 1):
            print(epoch, end='\r')
            # plt.cla()
            population = optimizer[epoch]
            population = optimizer.calc_rank(population)
            front = [x for x in population if x.rank == 1]

            indivs = [fit.get_indiv() for fit in population]

            if hasattr(population[0], 'id'):
                idx = max(*map(attrgetter('id'), population), idx)
            else:
                idx = max(*map(attrgetter('id'), indivs), idx)
            data.append([idx, len(front)])
        np.save(cache, data)


def ga_result2(out='result'):
    file = ut.fsort(glob.glob(os.path.join(out, f'epoch*.pickle')))[-1]
    optimizer = MOEAD.load(file=file)

    population = optimizer[-1]
    front = [x for x in population if x.rank == 1]
    front.sort(key=attrgetter('data.value'))

    for fit in front:
        # print(ind.value, ind.data.value)
        print(fit.get_indiv().get_gene().get_ivalue(),
              fit.value)

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

def __test__():
    print(ut.strnow('%Y%m%d_%H%M'))

def pltenv():
    rcParams['figure.figsize'] = 11, 7 # default=>(6.4, 4.8)
    rcParams['font.family'] = 'Times New Roman', 'serif'
    rcParams['font.size'] = 24 # default=>10.0

    # rcParams["mathtext.rm"] = 'Times New Roman'
    # rcParams["mathtext.it"] = 'Times New Roman'
    # rcParams["mathtext.bf"] = 'Times New Roman'
    # rcParams["mathtext.rm"] = 'Times New Roman'
    # rcParams["mathtext.sf"] = 'Times New Roman'
    # rcParams["mathtext.tt"] = 'Times New Roman'

    rcParams['savefig.directory'] = 'data'
    rcParams['savefig.transparent'] = True


def get_args():
    '''
    docstring for get_args.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('method', nargs='?', default='',
                        help='Main method type')
    parser.add_argument('--model', '-m', default='nsga2',
                        help='Model type')
    parser.add_argument('--out', '-o', default='',
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
    model = args.model
    out = os.path.join('result', model, args.out)
    clear = args.clear

    if args.test:
        __test__()
        return

    pltenv()

    if args.method == 'm1':
        ga_main1(model, out=out, clear_directory=clear)
    elif args.method == 'm2':
        ga_main2(out=out, clear_directory=clear)
    elif args.method == 'r1':
        ga_result1(out=out)
    elif args.method == 'r11':
        ga_result11(model, out=out)
    elif args.method == 'r11s':
        ga_result11(model, out=out, show=True)
    elif args.method == 'r12':
        ga_result12(model, out=out)
    elif args.method == 'r12s':
        ga_result12(model, out=out, show=True)
    elif args.method == 'r13':
        ga_result13(model, out=out)
    elif args.method == 'r13s':
        ga_result13(model, out=out, show=True)
    elif args.method == 'r14':
        ga_result14(model, out=out)
    elif args.method == 'r2':
        ga_result2(out=out)


if __name__ == '__main__':
    main()

'''
Note:
2018.11.21 08:22
MOEA/D ksize=50 N=500 # ksize大きくすると多様性が減少
2018.11.21 08:25
MOEA/D ksize=20 N=500 # ksize大きくすると多様性が減少

'''
