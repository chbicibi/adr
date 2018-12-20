import argparse
import os
import sys
import orbit_model as om

from matplotlib import rcParams
import utils as ut


def ga_main1(method_name, out='result', clear_directory=False):
    ''' GA実行
    '''
    if clear_directory and os.path.isdir(out):
        shutil.rmtree(out)

    popsize = 100
    epoch = 100
    ksize = 5
    save_trigger = lambda i: i == epoch # 最後だけ

    model_env = {'nsga2':om.NSGA2_ENV, 'moead':om.MOEAD_ENV}[method_name]

    with model_env(popsize=popsize, ksize=ksize) as env:
        optimizer = env.optimizer
        creator = env.creator
        with ut.stopwatch('main'):
            # GA開始
            # 初期集団生成
            population = optimizer.init_population(creator, popsize=popsize)
            history = [population]

            # 進化
            for i in range(1, epoch + 1):
                # optimizer.advance()
                population = optimizer(population)
                history.append(population)

                print('epoch:', i, 'popsize:', len(population), end='\r')
                if save_trigger(i):
                    file = f'popsize{popsize}_epoch{i}_{ut.strnow("%Y%m%d_%H%M")}.pkl'
                    file = os.path.join(out, file)
                    print('save:', file)
                    # optimizer.save(file=os.path.join(out, file))
                    ut.save(file, (env, optimizer, history))
            return env, optimizer, history


def get_plot():
    ###
    def ln(ind, origin):
        if not origin:
            return []
        return list(zip(origin[0].value, ind.value, origin[1].value))
    def dist(ind, origin):
        if not origin:
            return 0
        values = np.array([[ind.value, par.value] for par in origin])
        diff = values[:, :, 0] - values[:, :, 1]
        dists = np.sqrt(np.sum(diff ** 2, axis=1))
        return np.sum(dists)

    fig, ax = plt.subplots()
    def plot(pop):
        pairs = [(fit.data, fit.data.origin.origin or ()) for fit in pop]
        parents = list(chain(*(x[1] for x in pairs)))
        lines = [ln(ind, origin) for ind, origin in pairs]
        # print(lines)
        # exit()
        if False:
            if lines and [x for x in lines if x]:
                dists = [dist(ind, origin) for ind, origin in pairs]
                print(sum(dists))

        ax.cla()
        ax.set_xlim((0, 100))
        ax.set_ylim((0, 0.45))
        cm = plt.get_cmap('jet')
        # print(cm(10))
        # exit()

        if parents:
            x_p, y_p = np.array([ind.value for ind in parents]).T
            ax.scatter(x_p, y_p, c='r')

            for i, l in enumerate(lines):
                if l:
                    ax.plot(*l, c=cm(i/(len(lines)+1)), linewidth=0.5)

        x, y = np.array([fit.data.value for fit in pop]).T
        ax.scatter(x, y, c='b')
        # plt.pause(1e-10)
        ax.annotate("MOEA/D", xy=(0.5, -0.08), xycoords="axes fraction", fontsize=28, horizontalalignment="center", verticalalignment="top")
    plot.fig = fig
    return plot


def ga_res_temp1(out='result'):
    def get_model():
        # モデル読み込み
        # model_cls = {'nsga2':NSGA2, 'moead':MOEAD}[model]
        files = ut.fsort(glob.glob(os.path.join(out, f'*epoch*.pkl')))
        for i, file in enumerate(files):
            print(f'[{i}]', file)
        print('select file')
        n = int(input())
        if n < 0:
            return
        file = files[n]
        print('file:', file)
        env, optimizer, history = ut.load(file)
        return env, optimizer, history
        # optimizer = model_cls.load(file=file)

    def resume_main(env, optimizer, history):
        print('resume_main')
        plot = get_plot()
        for i, population in enumerate(history):
            plot(population)
            origin = population[0].data.origin.origin or []
            # print([x.id for x in origin], '->', population[0].data.id)
        plt.show()


    def resume_anim(env, optimizer, history):
        plot = get_plot()
        fig = plot.fig
        def update(i):
            print(i)
            population = history[i]
            # for i, population in enumerate(history):
            plot(population)
                # yield
        show_anim(fig, update, file='result/anim.mp4', frames=100)

    env, optimizer, history = get_model()
    print(type(history[0][0].get_indiv()))
    return
    resume_anim(env, optimizer, history)


def ga_res_temp2(out='result'):
    ''' クラス設計の確認用
    '''
    n_dim = 10 # 回収デブリの数
    Genome.set_weight([1, 3600] * n_dim)


    file = "result/nsga2/popsize100_epoch100_20181210_1318.pkl"
    env, optimizer, history = ut.load(file)
    pop0 = history[0] # 第0世代の解集団
    fit0 = pop0[0] # 解集団の0番目の個体の評価値
    print(type(fit0))

    ind0 = fit0.get_indiv() # 評価値の個体
    print(type(ind0))

    ''' 設計変数 '''
    genome = ind0.get_gene() # 個体が持つ遺伝子
    print(type(genome))

    dv_int = genome.get_ivalue() # デブリ番号
    dv_real = genome.get_dvalue() # 遷移時刻

    print(dv_int)
    print(dv_real)

    ''' 目的関数 '''
    dv, rcs = ind0.value # 評価値(生の値) => (delta_V [km/s], RCS [m^2])
    print(dv, rcs)


def ga_res_temp3(out='result'):
    ''' 結果プロットのテスト (散布図)
    '''
    n_dim = 10 # 回収デブリの数
    Genome.set_weight([1, 3600] * n_dim)
    file = "result/nsga2/popsize100_epoch100_20181210_1318.pkl"
    env, optimizer, history = ut.load(file)

    ''' プロット '''

    for pop in history:
        plt.cla()
        inds = [fit.get_indiv() for fit in pop]
        data = [ind.value for ind in inds] # [(dv0, rcs0), (dv1, rcs1), ...]
        data_np = np.array(data) # [[dv0, rcs0], [dv1, rcs1], ...]
        data_np_T = data_np.T # [[dv0, dv1, ...], [rcs0, rcs1, ...]]
        plt.scatter(*data_np_T)
        plt.pause(0.1)


def ga_res_temp4(out='result'):
    ''' 結果プロットのテスト (dvの履歴)
    '''
    n_dim = 10 # 回収デブリの数
    Genome.set_weight([1, 3600] * n_dim)
    file = "result/nsga2/popsize100_epoch100_20181210_1318.pkl"
    env, optimizer, history = ut.load(file)

    ''' プロット '''
    dvs = [[fit.get_indiv().value[0] for fit in pop] for pop in history] # [[dv0_0, dv0_1, ...], [dv1_0, dv1_1, ...], ...]
    dvs_min = [min(dv) for dv in dvs] # => [dvmin_0, dvmin_1, ...]

    plt.plot(dvs_min)
    plt.show()




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
    elif args.method == 'r1':
        ga_res_temp1(out=out)
    elif args.method == 'r2':
        ga_res_temp2(out=out)
    elif args.method == 'r3':
        ga_res_temp3(out=out)
    elif args.method == 'r4':
        ga_res_temp4(out=out)



if __name__ == '__main__':
    main()

'''
Note:
2018.11.21 08:22
MOEA/D ksize=50 N=500 # ksize大きくすると多様性が減少
2018.11.21 08:25
MOEA/D ksize=20 N=500 # ksize大きくすると多様性が減少

'''
