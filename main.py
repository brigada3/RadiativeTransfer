#!/usr/bin/env python3
# -*-coding: utf-8 -*-

import sys
import logging
import scipy.integrate
from b_system import b, Px, epsilon_x, kappa_x, F, f    


def complex_quadrature(func, a, b, **kwargs):
    def real_func(x):
        return scipy.real(func(x))
    def imag_func(x, *args):
        return scipy.imag(func(x))
    real_integral = scipy.integrate.quad(real_func, a, b, **kwargs)
    imag_integral = scipy.integrate.quad(imag_func, a, b, **kwargs)
    return (real_integral[0] + 1j*imag_integral[0], real_integral[1:], imag_integral[1:])


def main():
    sys.setrecursionlimit(10000000)
    logging.basicConfig(filename='log', filemode='w', level=logging.DEBUG)

    v, w0, mu, n, n1 = 1-0.5*1j, 0.95, 0.8, 20, 480

    for r in range(0, n+1):
        for i in range(0, n-r+1):
            b_value = b(alpha=i, v=v, w0=w0, mu_1=mu, n=n, n1=n1, r=r)
            logging.info('B[r={}][{}] = {}'.format(r, i, b_value))
           
        f_res = F(r=r, v=v, w0=w0, mu=mu, n=n, n1=n1)
        print (f_res)
        logging.info('F[r={}] = {}'.format(r, f_res))
        
        def math_func(fi):
            return f_res * scipy.exp(1j*fi*r)

        result = complex_quadrature(math_func, 0, 2*scipy.pi)[0]
        logging.info('Integrate F[r={}] = {}'.format(r, result))
   

if __name__ == '__main__':
    main()


 #if 1 <= i < n-r:
                #delta = (
                #    1j*v*(
                #        epsilon_x(i+1, r) * b(alpha=i+1, v=v, w0=w0, mu_1=mu, n=n, n1=n, r=r) +
                #        epsilon_x(i, r) * b(alpha=i-1, v=v, w0=w0, mu_1=mu, n=n, n1=n, r=r)
                #    ) - kappa_x(i, r, w0) * b(alpha=i, v=v, w0=w0, mu_1=mu, n=n, n1=n, r=r) + 
                #    (2*(r+i)+1)*f(r+i)*Px(r+i, r, mu)
                #).simplify()
                #logging.info('delta = {}'.format(delta))