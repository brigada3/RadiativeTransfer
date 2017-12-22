#!/usr/bin/env python3
# -*-coding: utf-8 -*-

import sys
import logging
from b_system import b, Px, epsilon_x, kappa_x, F, f    

def main():
    sys.setrecursionlimit(10000000)
    logging.basicConfig(filename='log', filemode='w', level=logging.INFO)

    v, w0, mu, r, n, n1 = 1-0.5j, 0.95, 0.8, 0, 20, 600

    for r in range(n):
        for i in range(n-r+1):
            logging.info('B[r={}][{}] = {}'.format(r, i, 
                b(alpha=i, v=v, w0=w0, mu_1=mu, n=n, n1=n1, r=r)))
            if 1 <= i < n-r:
                delta = (
                    1j*v*(
                        epsilon_x(i+1, r) * b(alpha=i+1, v=v, w0=w0, mu_1=mu, n=n, n1=n, r=r) +
                        epsilon_x(i, r) * b(alpha=i-1, v=v, w0=w0, mu_1=mu, n=n, n1=n, r=r)
                    ) - kappa_x(i, r, w0) * b(alpha=i, v=v, w0=w0, mu_1=mu, n=n, n1=n, r=r) + 
                    (2*(r+i)+1)*f(r+i)*Px(r+i, r, mu)
                ).simplify()
                logging.info('delta = {}'.format(delta))

    for i in range(n):
        logging.info('F[r={}] = {}'.format(i, F(i, v, w0, mu, n, n1)))
    
   

if __name__ == '__main__':
    main()
