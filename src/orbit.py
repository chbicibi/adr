#! /usr/bin/env python3

'''
Abstruct
'''

import argparse
import ctypes
import os
import shutil

import numpy as np


class Orbit(ctypes.Structure):
    _fields_ = [
        ('epc', ctypes.c_double),
        ('inc', ctypes.c_double),
        ('ran', ctypes.c_double),
        ('ecc', ctypes.c_double),
        ('ap', ctypes.c_double),
        ('ma', ctypes.c_double),
        ('mm', ctypes.c_double),
        ('type', ctypes.c_int32),
        ('sma', ctypes.c_double),
        ('slr', ctypes.c_double),
        ('prd', ctypes.c_double),
        ('axis_p', ctypes.c_double * 3),
        ('axis_q', ctypes.c_double * 3),
        ('axis_w', ctypes.c_double * 3),
        ('amom', ctypes.c_double * 3),
        ('dt', ctypes.c_double),
        ('ea', ctypes.c_double),
        ('ta', ctypes.c_double),
        ('pos', ctypes.c_double * 3),
        ('vel', ctypes.c_double * 3),
        ('dran', ctypes.c_double),
        ('dap', ctypes.c_double),
    ]

    def __init__(self):
        self.ptr = ctypes.byref(self)

    # def __init__(self):
    #     mod_name = 'orbit_base'
    #     proc_name = ''
    #
    #     cdll = np.ctypeslib.load_library('liborbit.dll', '.')
    #     f_add = getattr(cdll, f'__{mod_name}_MOD_{proc_name}')
    #     f_add.argtypes = [ctypes.POINTER(ctypes.c_double),
    #                       ctypes.POINTER(ctypes.c_double)]
    #     f_add.restype = ctypes.c_void_p


################################################################################


def __test__():
    point3d = ctypes.c_double*3
    print(point3d)
    print(type(point3d))
    print(dir(point3d))


def get_args():
    '''
    docstring for get_args.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', '-o', default='new_script',
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
        __test__()
        return

    file = args.out



if __name__ == '__main__':
    main()
