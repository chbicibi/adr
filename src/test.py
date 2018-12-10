# Python から Fortran のモジュールを呼ぶテスト

import ctypes
import numpy as np

def test0():
    '''
    Fortranサブルーチン呼び出し
    引数がスカラー値
    '''
    addmodule = ctypes.cdll.LoadLibrary('test.dll')

    # 引数の型を指定する。今回はintのポインタ型
    # Fortranのサブルーチン名"add"にアンスコ"_"をつける(gccでは)
    addmodule.add_.argtypes = [ctypes.POINTER(ctypes.c_int64),
                               ctypes.POINTER(ctypes.c_int64)]

    # 戻り値の型を指定する。Fortranのサブルーチンはvoidしか返せない
    addmodule.add_.restype = ctypes.c_void_p

    # 呼び出しに使う引数はctypesの型でラップする
    a, b = 10, 8
    a = ctypes.c_int64(a)
    b = ctypes.c_int64(b)

    # byrefでポインタにして渡す
    addmodule.add_(ctypes.byref(a), ctypes.byref(b))
    print(b.value) # 8

def test1():
    '''
    Fortranサブルーチン呼び出し
    引数がnumpy配列
    '''
    add_np = np.ctypeslib.load_library('test.dll', '.')
    add_np.add_array_.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        ctypes.POINTER(ctypes.c_int64)]
    add_np.add_array_.restype = ctypes.c_void_p

    a = np.arange(0, 10, dtype=np.float64) # 0,1,2,3,4,5,..,9
    b = a * 2 + 100                        # 0,2,4,6,8,...,18 (+100)
    size = ctypes.byref(ctypes.c_int64(b.size))
    add_np.add_array_(a+1, b, size)
    print(b) # [ 0 3 6 9 12 .. 27 ]

if __name__ == '__main__':
    test1()
