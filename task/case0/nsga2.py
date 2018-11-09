#! /usr/bin/env python3

'''
設計変数: 任意
目的関数: 1次元配列
解集団形状: 1次元配列, 2重配列(島モデル)

各個体にランク(自然数), 混雑度(実数)が割り当てられる
集団評価値(参考): ハイパーボリューム

解集団 =~ [個体]

選択
[個体] => [個体]

個体が生成されるのは次の3通りのいずれか
初期化: [] => 個体
交叉: [個体] => 個体 / [個体] => [個体]
突然変異: 個体 => 個体
'''

import argparse
import os
import pickle
import shutil

from functools import total_ordering

import numpy as np

from base import NondominatedSort


# class Initializer(object):
#     ''' [0, 1)の範囲の一様乱数による実数配列を返す '''

#     def __init__(self, size):
#         self._size = size

#     def __call__(self):
#         return np.random.rand(self._size)


@total_ordering
class Individual(object):
    ''' 進化計算個体
    '''
    current_id = 0

    def __init__(self, initial_gene=None):
        self.id = Individual.current_id
        if initial_gene is not None:
            self.gene = initial_gene
        else:
            self.gene = None
        self.value = None

        self.rank = None
        self.crowding = None

        Individual.current_id += 1

    def __call__(self):
        return self.value

    def __getitem__(self, key):
        return self.value[key]

    def __iter__(self):
        for v in self.value:
            yield v

    def __eq__(self, other):
        if not isinstance(other, Individual):
            return NotImplemented
        return self.fitness == other.fitness

    def __lt__(self, other):
        if not isinstance(other, Individual):
            return NotImplemented
        return self.fitness < other.fitness

    def __str__(self):
        return f'indiv{self.id}'

    def initialize(self, initializer):
        self.gene = initializer()

    def set_fitness(self, rank, fitness):
        self.rank = rank
        self.fitness = fitness

    def get_gene(self):
        ''' getter of genotype '''
        return self.gene

    def get_variable(self):
        ''' getter of phenotype '''
        return self.decode(self.gene)

    def encode(self, x):
        ''' phenotype -> genotype '''
        return x

    def decode(self, x):
        ''' genotype -> phenotype '''
        return x

    def evaluate(self, problem):
        if not self.value:
            self.problem = problem
            self.value = np.array(problem(self.get_variable()))


class Population(object):
    ''' 解集団
    GAの個体を管理する
    集団形状が1次元配列の場合はリストに類似
    島モデルの場合はmigrationを担当する
    '''

    def __init__(self, size, initializer=None):
        self.size = size
        self.initializer = initializer

        if initializer:
            self.values = [initializer() for i in range(size)]
        else:
            self.values = []

    def __call__(self):
        return self.values

    def __getitem__(self, key):
        return self.values[key]

    def __iter__(self):
        for v in self.values:
            yield v

    def __len__(self):
        return len(self.values)

    def append(self, x):
        self.values.append(x)

    def clear(self):
        self.values.clear()

    def filled(self):
        return len(self) >= self.size

################################################################################

class NSGA2(object):
    ''' NSGA-IIモデル '''

    def __init__(self, popsize, selection, crossover, mutation):
        self.popsize = popsize
        self.selection = selection
        self.crossover = crossover
        self.mutation = mutation

        self.initializer = None
        self.current_population = Population(popsize, Individual)
        self.next_population = Population(popsize)
        self.sort_fn = NondominatedSort(popsize)
        self.history = []

    def __call__(self):
        pass
        # self.advance(self.current_population)

    def __getitem__(self, key):
        return self.history[key]

    def __len__(self):
        return len(self.history)

    def set_initializer(self, initializer):
        ''' 初期集団生成用クラスを登録 '''
        self.initializer = initializer

    def setup(self, problem):
        ''' 最適化問題を登録 '''
        self.problem = problem

    def start(self):
        ''' 初期集団生成 '''
        for indiv in self.current_population:
            indiv.initialize(self.initializer)
            indiv.evaluate(self.problem)
        self.calc_fitness(self.current_population)
        self.history.append(self.current_population)

    def calc_fitness(self, population):
        ''' 各個体の集団内における適応度を計算する
        1. 比優越ソート
        '''
        rank = self.sort_fn(population)
        for ind, r in zip(population, rank):
            fit = 0.8 ** (r - 1)
            ind.set_fitness(r, fit)
        # for ind in self.current_population:
        #     print(ind.rank, ind.fitness)

    def advance(self):
        '''評価値計算→選択→交叉→突然変異→評価→世代交代'''
        # for ind in self.current_population:
        #     print(ind.fitness)

        # print('#1')

        while not self.next_population.filled():
            # print('#1.1')
            parents = self.get_parents(self.current_population)
            # print('#1.2')
            child = self.crossover(parents)
            # print('#1.3')
            child = self.mutation(child)
            # print('#1.4')
            child = Individual(child)
            # print('#1.5')
            child.evaluate(self.problem)
            self.next_population.append(child)

        self.calc_fitness(self.next_population)

        # print('#2')

        self.current_population = self.next_population
        self.history.append(self.current_population)
        self.next_population = Population(self.popsize)


    def get_parents(self, population):
        rest = list(population)
        res = []
        for i in range(2):
            parent, rest = self.selection(rest)
            res.append(parent)
        return res

    def get_elite(self):
        return [x for x in self.current_population if x.rank == 1]

    def save(self, file):
        dirname = os.path.dirname(file)
        if dirname:
            os.makedirs(dirname, exist_ok=True)
        with open(file, 'wb') as f:
            pickle.dump(self, f)

    def load(file):
        try:
            with open(file, 'rb') as f:
                return pickle.load(f)

        except FileNotFoundError:
            print('NSGA2.load: File is not found')
            raise


################################################################################

# def hypervolume


################################################################################

def __test__():
    initializer = Initializer(10)
    print(initializer())


def get_args():
    '''
    docstring for get_args.
    '''
    parser = argparse.ArgumentParser()
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

    file = args.out

    if not os.path.splitext(file)[1] == '.py':
        file = file + '.py'

    if args.force:
        pass
        # if os.path.exists(file):

    shutil.copy(__file__, file)
    print('create:', file)


if __name__ == '__main__':
    main()
