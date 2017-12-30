#!/usr/bin/env python3
# -*-coding: utf-8 -*-

import sys
import logging
import scipy
from b_system import F
from mpi import MpiProcess


w0, mu, n, n1 = 0.95, 0.8, 1000, 600


def calc_f_v_result(v):
    F_v_result = 0

    for r in range(0, n+1):
        F_r_result = F(r=r, v=v, w0=w0, mu=mu, n=n, n1=n1)

        fi, step = -scipy.pi, 0.1
        while -scipy.pi <= fi <= scipy.pi:
            F_v_result += F_r_result * scipy.exp(1j*r*fi)
            fi += step

        logging.info('F[r=%s][v=%s] = %s', r, v, F_r_result)

    return F_v_result


def main():
    sys.setrecursionlimit(10000000)

    F_v = {}

    def result_collector(mpi, data):
        v, i, F_v_result = data

        F_v[v] = F_v_result/(2*scipy.pi)

    def worker(mpi):
        while 0 == 0:
            data = mpi.recv_from_master()
            if data is None:
                break

            v, i = data

            F_v_result = calc_f_v_result(v)

            mpi.send_to_master((v, i, F_v_result))

    def master(mpi):
        v = -1 - 1j
        im_step, re_step = 0.1j, 0.5

        for i in range(101):
            v += im_step
            if not -1 < v.imag < 1:
                v = v - v.imag*1j + v.real + re_step

            mpi.send_to_worker((v, i))

        mpi.wait_all_workers()

        result = sum(F_v.values())
        logging.info('F_v sum: %s', result)

    process = MpiProcess(worker, master, result_collector)

    logging.basicConfig(
        filename='log' + str(process.rank),
        filemode='w',
        level=logging.INFO)

    process.run()


if __name__ == '__main__':
    main()


#for i in range(0, n-r+1):
    #   b_value = b(alpha=i, v=1-0.5*1j, w0=w0, mu_1=mu, n=n, n1=n1, r=r)
    #   logging.info('B[r={}][v={}][{}] = {}'.format(r, v, i, b_value))
