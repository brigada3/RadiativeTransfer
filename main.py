#!/usr/bin/env python3
# -*-coding: utf-8 -*-

import sys
import logging
import scipy
from b_system import b, F


def main():
    sys.setrecursionlimit(10000000)
    logging.basicConfig(filename='log', filemode='w', level=logging.DEBUG)

    w0, mu, n, n1 = 0.95, 0.8, 1000, 600

    F_v = {}
    v = -1 - 1j
    im_step, re_step = 0.1j, 0.5
    for _ in range(101):
        v += im_step
        if not -1 < v.imag < 1:
            v = v - v.imag*1j + v.real + re_step

        F_v_result = 0
        for r in range(619, n+1):
            F_r_result = F(r=r, v=v, w0=w0, mu=mu, n=n, n1=n1)

            fi, step = -scipy.pi, 0.1
            while -scipy.pi <= fi <= scipy.pi:
                F_v_result += F_r_result * scipy.exp(1j*r*fi)
                fi += step

            logging.info('F[r=%s][v=%s] = %s', r, v, F_r_result)

        F_v[v] = F_v_result/F_v_result
        logging.info('F[v=%s] = %s', v, F_v_result)

    result = sum(F_v.values)
    logging.info('F_v sum: %s', F_v_result)



if __name__ == '__main__':
    main()


#for i in range(0, n-r+1):
    #   b_value = b(alpha=i, v=1-0.5*1j, w0=w0, mu_1=mu, n=n, n1=n1, r=r)
    #   logging.info('B[r={}][v={}][{}] = {}'.format(r, v, i, b_value))
