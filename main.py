#!/usr/bin/env python3
# -*-coding: utf-8 -*-

import sys
import logging
from b_system import b, Px, epsilon_x, kappa_x, F, f    

def main():
    sys.setrecursionlimit(10000000)

    v, w0, mu, r, n, n1 = 1-0.5j, 0.95, 0.8, 0, 20, 600

    b_values = []
    for i in range(21):
        b_values.append(b(i, v, w0, mu, r, n, n1))
        print('b[{}] = {}'.format(i, b_values[i]))
        print('F[{}] = {}'.format(i, F(i, r, v, w0, mu, n, n1)))

    for i in range(0, n-r):
        delta = (
            1j*v*(
                epsilon_x(i+1, r)*b_values[i+1] +
                epsilon_x(i, r)*b_values[i-1]
            ) - kappa_x(i, r, w0) * b_values[i] + (2*(r+i)+1)*f(r+i)*Px(r+i, r, mu)
        ).simplify()
        print('delta', delta)


if __name__ == '__main__':
    main()
