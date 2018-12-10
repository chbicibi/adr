import pickle


class ClassName0(object):
    ''' デフォルト '''
    def __init__(self, arg):
        print('init:', arg)
        self.arg = arg
        self.message = 'arg=' + arg


class ClassName1(object):
    ''' スタンダード '''
    def __new__(cls, arg):
        print('new:', arg)
        return super().__new__(cls)

    def __init__(self, arg):
        print('init:', arg)
        self.arg = arg
        self.message = 'arg=' + arg

    def __reduce_ex__(self, protocol):
        return type(self), (self.arg,)


class ClassName2(object):
    ''' __init__が呼ばれない
    __new__がある場合は__getnewargs__は必須
    '''
    def __new__(cls, arg):
        print('new:', arg)
        return super().__new__(cls)

    def __init__(self, arg):
        print('init:', arg)
        self.arg = arg
        self.message = 'arg=' + arg

    def __getnewargs__(self):
        return self.arg,


class ClassName3(object):
    ''' __setstate__を定義した場合は__getstate__の返す値は自由
    '''
    # def __new__(cls, arg):
    #     print('new:', arg)
    #     return super().__new__(cls)

    def __init__(self, arg):
        print('init:', arg)
        self.arg = arg
        self.message = 'arg=' + arg

    def __getstate__(self):
        print('__getstate__')
        return {'arg': self.arg, 'message': self.message}

    def __setstate__(self, state):
        print('__setstate__')
        self.__dict__ = state


def __test__(ClassName):
    obj = ClassName('hello')


def __test_save__(ClassName):
    obj = ClassName('hello')
    with open('test.pkl', 'wb') as f:
        pickle.dump(obj, f)


def __test_load__(ClassName):
    with open('test.pkl', 'rb') as f:
        obj = pickle.load(f)
        print(obj.message)


def main():
    print('#0')
    __test__(ClassName0)
    print('#1')
    __test_save__(ClassName0)
    print('#2')
    __test_load__(ClassName0)


if __name__ == '__main__':
    main()
