#!/usr/bin/env python3
# -*-coding: utf-8 -*-

import sys
import sympy
import b_system

def main():
    sympy.init_printing()
    sys.setrecursionlimit(100000)

    s, r, w0, i, n = sympy.symbols('s r w0 i n')
    alpha, l, j, v = sympy.symbols('alpha l j v')
    mu, mu_1 = sympy.symbols('mu mu_1')

    for i in range(0, 21):
        b_system.b(sympy.Integer(str(i)), -0.5j, 0.95, 0.8, sympy.Integer('0'))

if __name__ == '__main__':
    main()
