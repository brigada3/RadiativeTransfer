#!/usr/bin/env python3
# -*-coding: utf-8 -*-

import sympy, scipy.special
import functools
import logging
#logging.basicConfig(filename='log', filemode='w', level=logging.DEBUG)

G = 0.999


@functools.lru_cache(maxsize=None)
def f(s, g=G):
    return g**s


@functools.lru_cache(maxsize=None)
def epsilon_x(s, r):
    logging.debug('Enter epsilon_x')
    logging.debug('Epsilon_x args s={} r={}'.format(s, r))
    
    math_func = (s*(s+2*r))**(1/2)
    
    logging.debug('Epsilon_x args s={} r={}'.format(s, r))
    logging.debug('Epsilon_x result: {}'.format(math_func))    
    return math_func


@functools.lru_cache(maxsize=None)
def kappa_x(s, r, w0):
    logging.debug('Enter kappa_x')
    logging.debug('Kappa_x args: s={} r={} w0={}'.format(s,r,w0))
    
    math_func = (2*(s+r)+1) * (1-w0*f(s+r))

    logging.debug('Kappa_x args: s={} r={} w0={}'.format(s,r,w0))
    logging.debug('Kappa_x result {}'.format(math_func))
    return math_func


@functools.lru_cache(maxsize=None)
def fi(v_arg):
    logging.debug('Enter fi')
    logging.debug('fi args: v={}'.format(v))
    
    v = sympy.Symbol('v')
    math_func = sympy.Piecewise(
        (sympy.atan(sympy.im(v-1j)/sympy.re(v-1j)), sympy.re(v-1j) > 0),
        (-sympy.pi/2, sympy.And(sympy.Eq(sympy.re(v-1j), 0), sympy.im(v-1j) < 0)),
        (-sympy.pi/2 - sympy.atan(sympy.re(v-1j)/sympy.im(v-1j)), sympy.And(sympy.re(v-1j)<0, sympy.im(v-1j)<0)),
        (-sympy.pi, sympy.And(sympy.re(v-1j)<0, sympy.Eq(sympy.im(v-1j), 0))),
        (-sympy.pi + sympy.atan(sympy.im(v-1j)/sympy.re(v-1j)), sympy.And(sympy.re(v-1j)<0, sympy.im(v-1j)>0))

    )

    result = math_func.subs(v, v_arg)
    logging.debug('fi args: v={}'.format(v))
    logging.debug('fi result {}'.format(result))    
    return result


@functools.lru_cache(maxsize=None)
def fi_1(v_arg):
    v = sympy.Symbol('v')    
    math_func = sympy.Piecewise(
        (sympy.atan(sympy.im(v+1j)/sympy.re(v+1j)), sympy.re(v+1j) > 0),
        (sympy.pi/2, sympy.And(sympy.Eq(sympy.re(v+1j), 0), sympy.im(v+1j) > 0)),
        (sympy.pi/2 - sympy.atan(sympy.re(v+1j)/sympy.im(v+1j)), sympy.And(sympy.re(v+1j)<0, sympy.im(v+1j)>0)),
        (sympy.pi, sympy.And(sympy.re(v+1j)<0, sympy.Eq(sympy.im(v+1j), 0))),
        (sympy.pi + sympy.atan(sympy.im(v+1j)/sympy.re(v+1j)), sympy.And(sympy.re(v+1j)<0, sympy.im(v+1j)<0))
    )

    return math_func.subs(v, v_arg)


@functools.lru_cache(maxsize=None)
def q_x(s, r, w0, n):
    logging.debug("Enter q_x, called by gamma_x")
    logging.debug("Args: s={} r={} w0={}".format(s, r, w0))
    if s <= n-r-1:
        result = ((epsilon_x(s+1, r))**2) / (
            (2*(s+r+1)+1) * (1-w0*f(s+r+1)) * 
            (2*(s+r)+1) * (1-w0*f(s+r))
        )
    elif s == n-r:
        result = ((epsilon_x(s+1, r))**2) / (
            (2*(s+r)+1)*(1-w0*f(s+r)) *
            (2*(s+r+1)+1)
        )
    elif s >= n-r+1:
        result = ((epsilon_x(s+1, r))**2) / (
            (2*(s+r+1)+1) * (2*(s+r)+1)
        )
   
    logging.debug("Args: s={} r={} w0={}".format(s, r, w0))
    logging.debug("q_x result: {}".format(result))
    return result


@functools.lru_cache(maxsize=None)
def maple_f(i, v, r, n, n1):
    logging.debug('Enter maple_f')
    logging.debug('Args i={}, v={}, r={}, n={}, n1={}'.format(i, v, r, n, n1))
    
    if i == 0:
        A = (0.5)*(-1 + (1+v)**(0.5))
        result = (
            1 + a(l=0, v=v, r=r, n=n, n1=n1)/(1+A)
        )
    else:
        result = (
            1 + a(l=i, v=v, r=r, n=n, n1=n1)/maple_f(i-1, v,r,n,n1)
        )

    logging.debug('Args i={}, v={}, r={}, n={}, n1={}'.format(i, v, r, n, n1))
    logging.debug('maple_f result: {}'.format(result))
    return result


@functools.lru_cache(maxsize=None)
def gamma(l, v, r, n, n1):
    logging.debug('Enter gamma')
    logging.debug('Gamma args: l={}, v={}, r={}, n={}, n1={}'.format(l, v, r,  n, n1))

    result = maple_f(n1-1, v, r, n, n1)
    
    logging.debug('Gamma args: l={}, v={}, r={}, n={}, n1={}'.format(l, v, r, n, n1))    
    logging.debug("Gamma result: {}".format(result))
    return result


@functools.lru_cache(maxsize=None)
def gamma_x(s, v, r, w0, n, n1):
    logging.debug('Enter gamma_x')
    logging.debug('Gamma_x args: s={}, v={}, r={}, w0={}, n={}, n1={}'.format(s, v, r, w0, n, n1))
    
    if s == n-r+1:
        result = gamma(l=1, v=v, r=r, n=n, n1=n1)
        logging.debug('Gamma_x args: s={}, v={}, r={}, w0={}, n={}, n1={}'.format(s, v, r, w0, n, n1))        
        logging.debug("Gamma_x s==n-r+1 result: {}".format(result))
        return result
    elif isinstance(s, int):
        gamma_x_s_plus_1 = gamma_x(s+1, v, r, w0, n, n1)
    else:
        raise ValueError('s must be a number!')

    result = 1 + q_x(s, r, w0, n)*v/gamma_x_s_plus_1
    
    logging.debug('Gamma_x args: s={}, v={}, r={}, w0={}, n={}, n1={}'.format(s, v, r, w0, n, n1))    
    logging.debug("Gamma_x result: {}".format(result))
    return result


@functools.lru_cache(maxsize=None)
def h_x(alpha, v, r, w0, n, n1):
    logging.debug("Enter h_x")
    logging.debug("Args h_x alpha={} v={} r={} w0={}".format(alpha, v, r, w0))
    
    math_func = (
        1j*v*epsilon_x(s=alpha, r=r) /
        (kappa_x(s=alpha, r=r, w0=w0) * gamma_x(s=alpha, v=v**2, r=r, w0=w0, n=n, n1=n1))
    )

    logging.debug("Args h_x alpha={} v={} r={} w0={}".format(alpha, v, r, w0))
    logging.debug("H_x return: {}".format(math_func))
    return math_func


@functools.lru_cache(maxsize=None)
def a(l, v, r, n, n1):
    logging.debug('Enter a')
    logging.debug('a args: l={}, v={}, r={}, n={}'.format(l, v, r,  n))
    
    math_func = ((n+n1-l)**2 - r**2)*v/( (2*(n+n1-l)+1) * (2*(n+n1-l+1)+1) )
    
    logging.debug('a args: l={}, v={}, r={}, n={}'.format(l, v, r,  n))
    logging.debug('a result: {}'.format(math_func))
    return math_func


@functools.lru_cache(maxsize=None)
def Px(i, r, mu):
    logging.debug('Enter Px')
    logging.debug('Px args: i={}, r={}, mu={}'.format(i, r, mu))
    if i == r:
        result = (
            (1-mu**2)**(0.5*r) * 2**(-r) *
            (scipy.special.gamma(2*r+1)**(0.5)) / scipy.special.gamma(r+1)
        )
    elif i == r + 1:
        result = (
            mu*(1-mu**2)**(0.5*r) * 2**(-r) *
            scipy.special.gamma(2*r+1)**(0.5) * (2*r+1)**(0.5) / scipy.special.gamma(r+1)
        )
    elif i >= 2 + r:
        result = (
            (mu*(2*i-1)*Px(i-1, r, mu) - 
            Px(i-2, r, mu) * ((i-1)**2 - r**2)**(0.5)) / 
            (i**2-r**2)**(0.5)
        )
    
    logging.debug('Px args: i={}, r={}, mu={}'.format(i, r, mu))    
    logging.debug('Px result: {}'.format(result))
    print 
    return float(result)


@functools.lru_cache(maxsize=None)
def psi_x(j, v, r, w0, n, n1):
    logging.debug('Enter psi_x called by b')
    logging.debug('Args psi_x j={} v={} r={} w0={} n={} n1={}'.format(j, v, r, w0, n, n1))
    
    result = 1
    if j != 0:
        for l in range(1, j+1):
            result *= h_x(alpha=l, v=v, r=r, w0=w0, n=n, n1=n1)

    logging.debug('Args psi_x j={} v={} r={} w0={} n={} n1={}'.format(j, v, r, w0, n, n1))
    logging.debug('psi_x result: {}'.format(result))
    return result


@functools.lru_cache(maxsize=None)
def nu_x(alpha, v, r, w0, mu_1, n, n1):
    logging.debug('Enter nu_x')
    logging.debug('Args alpha=%s v=%s r=%s w0=%s mu_1=%s n=%s n1=%s' % (alpha, v, r, w0, mu_1, n, n1))

    if alpha == n-r+1:
        logging.debug('Args alpha=%s v=%s r=%s w0=%s mu_1=%s n=%s n1=%s' % (alpha, v, r, w0, mu_1, n, n1))
        result = 0
        logging.debug('Nu_x result: %s' % result)        
        return result
    elif alpha == n-r:
        logging.debug('Condition alpha == n-r')
        
        result = (
            (2*n+1)*f(n)*Px(n, r, mu_1) /
            (kappa_x(n-r, r, w0)*gamma_x(s=n-r, v=v**2, r=r, w0=w0, n=n, n1=n1))
        )
        logging.debug('Args alpha=%s v=%s r=%s w0=%s mu_1=%s n=%s n1=%s' % (alpha, v, r, w0, mu_1, n, n1))
        logging.debug('Nu_x result: %s' % result)
        return result
    elif isinstance(alpha, int):
        logging.debug('Nu_x recirsuon alpha+=1')
        nu_x_alpha_plus_1 = nu_x(alpha+1, v, r, w0, mu_1, n, n1)
    else:
        ValueError('alpha')

    result = (
        (kappa_x(s=alpha, r=r, w0=w0)*gamma_x(s=alpha, v=v**2, r=r, w0=w0, n=n, n1=n1))**(-1) *
        (1j*v*epsilon_x(s=alpha+1, r=r)*nu_x_alpha_plus_1 + 
            (2*(alpha+r)+1)*f(alpha+r)*Px(alpha+r, r, mu_1)
        )
    )

    logging.debug('Args alpha=%s v=%s r=%s w0=%s mu_1=%s n=%s n1=%s' % (alpha, v, r, w0, mu_1, n, n1))
    logging.debug('Nu_x result: %s' % result)
    return result


@functools.lru_cache(maxsize=None)
def b(alpha, v, w0, mu_1, r, n, n1):
    logging.debug('Enter b')
    logging.debug('Args b: alpha={} v={} w0={} mu_1={} r={} n={} n1={}'.format(alpha, v, w0, mu_1, r, n, n1))
    
    if alpha == 0 and n-r >= 0:
        sum_result = 0
        for j in range(0, n-r+1):
            sum_result += psi_x(j=j, v=v, r=r, w0=w0, n=n, n1=n1)*(2*(j+r)+1)*f(j+r)*Px(j+r, r, mu_1)
        logging.debug('sum_result = {}'.format(sum_result))

        result = (
            (kappa_x(0, r, w0) * 
            gamma_x(0, v**2, r, w0, n, n1))**(-1) *
            sum_result
        )
        logging.debug('Args b: alpha={} v={} w0={} mu_1={} r={} n={} n1={}'.format(alpha, v, w0, mu_1, r, n, n1))
        logging.debug('b_0 result: {}'.format(result))
        return result
    elif isinstance(alpha, int) and 1 <= alpha <= n-r and n-r >= 1:
        b_alpha_minus_1 = b(alpha-1, v, w0, mu_1, r, n, n1)
    else:
        print(alpha)
        raise ValueError('alpha nust be integer')

    result = (
        h_x(alpha, v, r, w0, n, n1) * b_alpha_minus_1 +
        nu_x(alpha, v, r, w0, mu_1, n, n1)
    )
    
    logging.debug('Args b: alpha={} v={} w0={} mu_1={} r={} n={} n1={}'.format(alpha, v, w0, mu_1, r, n, n1))
    logging.debug('B result: %s' % result)
    return result


def F(r, v, w0, mu, n, n1):
    logging.debug('Enter F')
    logging.debug('Args F r={} v={} w0={} n={} n1={}'.format(r, v, w0, n, n1))
    
    result = 0
    for k in range(0, n-r+1):
        result += (
            (2*k+2*r+1) *
            b(alpha=k, v=v, w0=w0, mu_1=mu, r=r, n=n, n1=n1) *
            Px(k+r, r, mu)
        )
    result = result/2

    logging.debug('Args F r={} v={} w0={} n={} n1={}'.format(r, v, w0, n, n1))
    logging.debug('F result: {}'.format(result))
    return result
