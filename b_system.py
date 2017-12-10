#!/usr/bin/env python3
# -*-coding: utf-8 -*-

import sympy
import functools
import logging
#logging.basicConfig(filename='log', filemode='w', level=logging.DEBUG)

N = 20
N1 = 24*N
G = 0.999
B = []

def f(s, g=G):
    return g**s


def epsilon_x(s, r):
    return sympy.sqrt(s*(s+2*r))


def kappa_x(s, r, w0, caller):
    logging.debug('Enter kappa_x called by {}'.format(caller))
    logging.debug('Kappa_x args: s={} r={} w0={}'.format(s,r,w0))    
    math_func = (2*(s+r)+1) * (1-w0*f(s+r))
    logging.debug('Kappa_x return {}'.format(math_func))    
    return math_func


def fi(v_arg):
    v = sympy.Symbol('v')
    math_func = sympy.Piecewise(
        (sympy.atan(sympy.im(v-1j)/sympy.re(v-1j)), sympy.re(v-1j) > 0),
        (-sympy.pi/2, sympy.And(sympy.Eq(sympy.re(v-1j), 0), sympy.im(v-1j) < 0)),
        (-sympy.pi/2 - sympy.atan(sympy.re(v-1j)/sympy.im(v-1j)), sympy.And(sympy.re(v-1j)<0, sympy.im(v-1j)<0)),
        (-sympy.pi, sympy.And(sympy.re(v-1j)<0, sympy.Eq(sympy.im(v-1j), 0))),
        (-sympy.pi + sympy.atan(sympy.im(v-1j)/sympy.re(v-1j)), sympy.And(sympy.re(v-1j)<0, sympy.im(v-1j)>0))
    )
    return math_func.subs(v, v_arg)


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
def q_x(s_arg, r_arg, w0_arg, n_arg=N):
    #logging.debug("Enter q_x, called by gamma_x")    
    #logging.debug("Args: ", s_arg, r_arg, w0_arg)

    s, r, n, w0 = sympy.symbols('s r n w0')
    math_func = sympy.Piecewise(
        ((s*(s+2*r))/((2*(s+r)+1)*(2*(s+r+1)+1)*(1-w0*f(s+r)*(1-w0*f(s+r+1)))), sympy.And(s+r+1<=n, n>=2)),
        ((s*(s+2*r))/((2*(s+r)+1)*(2*(s+r+1)+1)*(1-w0*f(s+r))), sympy.Eq(s+r, n)),
        ((s*(s+2*r))/((2*(s+r)+1)*(2*(s+r+1)+1)), s >= n-r+1)
    )
    #logging.debug("Result q_x: ", math_func.subs({s:s_arg, r:r_arg, n:n_arg, w0:w0_arg}), end='\n\n')

    return math_func.subs({s:s_arg, r:r_arg, n:n_arg, w0:w0_arg})


@functools.lru_cache(maxsize=None)
def gamma(l, v, r, n=N, n1=N1):
    logging.debug('Enter gamma')
    logging.debug('Gamma args: l={}, v={}, r={}, n={}, n1={}'.format(l, v, r,  n, n1))
    
    if l == n1 + 1:
        result = (
            (-1 +
                (sympy.exp(1)**(fi(v) + fi_1(v))*1j/2) *
                ( (1+sympy.re(v)**2-sympy.im(v)**2)**2 + 4*sympy.re(v)**2*sympy.im(v)**2)**0.5
            ) / 2 + 1
        )
        logging.debug("Gamma result = {}".format(result))
        return result
    
    result = 1 + a(l, v, r, n)/gamma(l+1, v, r, n, n1)

    logging.debug("Gamma result: {}".format(result))
    return result


#def gamma_x(s, v, r, w0, n=N, n1=N1):
#    logging.debug("Enter gamma_x")
#    logging.debug("Args gamma_x: s={} r={} v={}".format(s, r, v))    
#    result = gamma(1, v, r, n, n1)
#    logging.debug("Init result via gamma = {}".format(result))
#
#    logging.debug("start ind gamma_x iteration")
#    for ind in range(n-r, s-1, -1):
#        #logging.debug("ind: ", ind)
#        result = 1 + q_x(ind, r, w0, n)*v/result
#        #logging.debug("new reuslt: ", result)
#    
#    logging.debug("Final result: {}".format(result))
#    
#    return result


@functools.lru_cache(maxsize=None)
def gamma_x(s, v, r, w0, n=N, n1=N1):
    if s == n-r+1:
        return gamma(1, v, r, n, n1)
    elif isinstance(s, int):
        gamma_x_s_plus_1 = gamma_x(s+1, v, r, w0, n, n1)
    else:
        raise ValueError('s must be a number!')

    return 1 + q_x(s, r, w0, n)*v/gamma_x_s_plus_1


@functools.lru_cache(maxsize=None)
def h_x(alpha, v, r, w0, caller, n=N, n1=N1):
    logging.debug("Enter h_x called by {}".format(caller))
    logging.debug("Args h_x alpha={} v={} r={} w0={}".format(alpha, v, r, w0))
    math_func = (
        1j*v*epsilon_x(alpha, r) *
        (kappa_x(alpha, r, w0, "h_x") * gamma_x((alpha), v**2, r, w0, n, n1))**(-1)
    )
    logging.debug("H_x return: {}".format(math_func))    
    return math_func


def a(l, v, r, n=N):
    math_func = ((n+l)**2-r**2)*v**2 /( (2*(n+l)+1) * (2*(n+l+1)+1) )
    return math_func


@functools.lru_cache(maxsize=None)
def Px(i, r, mu):
    logging.debug('Enter Px')
    logging.debug('Px args: i={}, r={}, mu={}'.format(i, r, mu))    
    
    if i == r-1:
        return 0
    elif i == r:
        return (1-mu**2)**(r/2) * (2)**(-r) * sympy.sqrt(sympy.factorial(2*r)/sympy.factorial(r)**2)
    elif isinstance(i, int):
        Px_i_minus_1 = Px(i-1, r, mu)
        Px_i_minus_2 = Px(i-2, r, mu)
    else:
        raise ValueError('i must be integer')
        
    result = (
        ( 2*i*mu*Px_i_minus_1 - sympy.sqrt((i-r-1)*(i+r-1)) * Px_i_minus_2 ) *
        ( (i-r)*(i+r) )**(-1/2)
    )
    logging.debug('Px result: {}'.format(result))

    return result


@functools.lru_cache(maxsize=None)
def psi_x(j, v, r, w0, n=N, n1=N1):
    logging.debug('Enter psi_x called by b')
    logging.debug('Args psi_x j={} v={} r={} w0={} n={} n1={}'.format(j, v, r, w0, n, n1))
    
    result = 1
    if j != 0:
        for l in range(1, j+1):
            result *= h_x(l, v, r, w0, "psi_x", n, n1)

    if not isinstance(result, int):
        result = result.simplify()

    logging.debug('psi_x return: {}'.format(result))
    return result


@functools.lru_cache(maxsize=None)
def nu_x(alpha, v, r, w0, mu_1, n=N, n1=N1):
    if alpha == n-r+1:
        return 0
    elif alpha == n-r:
        return (
            ((2*n+1)*f(n)*Px(alpha+r, r, mu_1)) /
            (kappa_x(n-r, r, w0, "nu_x")*gamma_x(n-r, v**2, r, w0, n, n1))
        )
    elif isinstance(alpha, int):    
        nu_x_alpha_plus_1 = nu_x(alpha+1, v, w0, mu_1, n, n1)
    else:
        ValueError('alpha')

    return (
        (kappa_x(alpha, r, w0, "nu_x")*gamma_x(alpha, v**2, r, w0, n, n1))**(-1)*
        (1j*v*epsilon_x(alpha+1, r)*nu_x_alpha_plus_1 + 
            (2*(alpha+r)+1)*f(alpha+r)*Px(alpha+r, r, mu_1)
        )
    )

@functools.lru_cache(maxsize=None)
def b(alpha, v, w0, mu_1, r, n=N, n1=N1):
    logging.debug('Enter b')
    logging.debug('Args b: alpha={} v={} w0={} mu_1={} r={} n={} n1={}'.format(alpha, v, w0, mu_1, r, n, n1))
    if alpha == 0 and n-r >= 0:
        sum_result = 0
        for j in range(0, n-r+1):
            sum_result += psi_x(j, v, r, w0, n, n1)*(2*(j+r)+1)*f(j+r)*Px(j+r, r, mu_1)
        logging.debug('sum_result = {}'.format(sum_result))  
        result = (
            (kappa_x(0, r, w0, "b") * 
            gamma_x(0, v**2, r, w0, n, n1))**(-1) *
            sum_result
        )
        logging.debug('b_0 = {}'.format(result))
        B.insert(0, result)
        return result
    elif isinstance(alpha, int) and 1 <= alpha <= n-r and n-r > 1:
        if alpha-1 < len(B):
            b_alpha_minus_1 = B[alpha-1]
        else:
            b_alpha_minus_1 = b(alpha-1, v, w0, mu_1, r, n, n1)
    else:
        raise ValueError('alpha nust be integer')

    result = (
        h_x(alpha, v, r, w0, "b", n, n1) * b_alpha_minus_1 +
        nu_x(alpha, v, r, w0, mu_1, n, n1)
    )
    B.insert(alpha, result)
    return result


def F(s, r, v, w0, mu, mu_1,  n=N, n1=N1):
    math_func = (0.5)*sympy.Sum(
        (2*(s+r+1))*b(s, v, w0, mu_1, r, n, n1)*Px(s+r, r, mu),
        (s, 0, n)
    )

    return math_func
