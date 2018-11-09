#! /usr/bin/env python3

'''
Abstruct
'''

import argparse
import os
import shutil


class Problem(object):
    ''' 最適化対象の問題 '''
    def __init__(self):
        pass

    def __call__(self, x):
        return x


################################################################################


def __test__():
    pass


def get_args():
    '''
    docstring for get_args.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', '-o', default='new_script.py',
                        help='Filename of the new script')
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
        __path__()
        return

    file = args.out

    if not os.path.exists(file):
        shutil.copy(__file__, file)
        print('create:', file)


if __name__ == '__main__':
    main()
