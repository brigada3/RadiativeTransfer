#!/usr/bin/env python3
# -*-coding: utf-8 -*-

import sys
import sympy
import b_system

def main():
    sympy.init_printing()

    s, r, w0, i, n = sympy.symbols('s r w0 i n')
    alpha, l, j, v = sympy.symbols('alpha l j v')
    mu, mu_1 = sympy.symbols('mu mu_1')

    sympy.pprint(b_system.gumma(sympy.Integer('480240'), 1+2j, r))
    sympy.pprint(b_system.gamma(sympy.Integer('480240'), 1+2j, r))
    

if __name__ == '__main__':
    main()
