#!/usr/bin/env python3
# -*-coding: utf-8 -*-

from sympy import *

def main():
    init_printing()



class h(Function):
    @classmethod
    def eval(cls, n):
        if n == 0:
            return 1
        else:
            return diff(h(x-1))


if __name__ == '__main__':
    main()
