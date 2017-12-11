#!/usr/bin/env python3
# -*-coding: utf-8 -*-

import sys
import logging
import sympy
import b_system

def main():
    sympy.init_printing()
    sys.setrecursionlimit(10000000)

    s, r, w0, i, n = sympy.symbols('s r w0 i n')
    alpha, l, j, v = sympy.symbols('alpha l j v')
    mu, mu_1 = sympy.symbols('mu mu_1')

    for i in range(21):
        print('b[{}] = {}'.format(i, b_system.b(i, 1-0.5j, 0.95, 0.8, 0)))


if __name__ == '__main__':
    main()
