#!/usr/bin/env python3
# -*-coding: utf-8 -*-

import sys
import logging
import b_system

def main():
    sys.setrecursionlimit(10000000)

    v, w0, mu_1, r = 1-0.5j, 0.95, 0.8, 0
    b = []

    for i in range(21):
        b.append(b_system.b(i, v, w0, mu_1, r))
        print('b[{}] = {}'.format(i, b[i]))
    
    for s in range(20): 
        left = (i*v*b_system.epsilon_x(s+1, r)*b[s+1] + i*v*b_system.epsilon_x(s-1, r)*b[s-1]).simplify()
        right = (b_system.kappa_x(s, r, w0)*b[s] - 
            (2*(s+r)+1)*b_system.f(s+r) * 
            b_system.Px(s+r, r, mu_1)).simplify()
        print('delta', left-right)


if __name__ == '__main__':
    main()
