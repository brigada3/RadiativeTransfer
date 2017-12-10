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

    print('b[2] = ', b_system.b(2, -0.5j, 0.95, 0.8, 0))

if __name__ == '__main__':
    main()
