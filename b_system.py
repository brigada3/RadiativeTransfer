#!/usr/bin/env python3
# -*-coding: utf-8 -*-

"""
    N   -- natural numbers
    N0  -- natural numbers w0ith 0
"""


import math


F = []
N = len(F)
# TODO: Find out what is n1
N1 = 42


def epsilon_x(s, r):
    r = abs(r)
    return (s*(s+2*r))**0.5


def kappa_x(s, r, w0, f=F):
    r = abs(r)
    return(2*(s+r)+1)*(1-w0*f[s+r])


def fi(v):
    v = v - 1j

    if v.real > 0:
        return math.atan(v.imag/v.real)
    elif v.real == 0 and v.imag < 0:
        return -math.pi/2
    elif v.real < 0 and v.imag < 0:
        return -math.pi/2 - math.atan(v.real/v.imag)
    elif v.real < 0 and v.imag == 0:
        return -math.pi
    elif v.real < 0 and v.imag > 0:
        return -math.pi + math.atan(v.imag/v.real)


def fi_1(v):
    v = v + 1j

    if v.real > 0:
        return math.atan(v.imag/v.real)
    elif v.real == 0 and v.imag < 0:
        return math.pi/2
    elif v.real < 0 and v.imag < 0:
        return math.pi/2 - math.atan(v.real/v.imag)
    elif v.real < 0 and v.imag == 0:
        return math.pi
    elif v.real < 0 and v.imag > 0:
        return math.pi + math.atan(v.imag/v.real)


def q_x(s, r, w0, n=N, f=F):
    r = abs(r)
    result = s*(s+2*r)

    if s+1+r <= n and n >= 2:
        return result/((2*(s+r)+1)*(2*(s+r+1)+1)*(1-w0*f[s+r])*(1-w0*f[s+r+1]))
    elif s+r == n:
        return result/((2*(s+r)+1)*(2*(s+r+1)+1)*(1-w0*f[s+r]))
    elif s >= n-r+1:
        return result/((2*(s+r)+1)*(2*(s+r+1)+1))


def h_x(alpha, v, r, w0):
    """
        alpha belong N
        asb(r) belong N0
    """
    r = abs(r)
    return 1j*v*epsilon_x(alpha, r)/(kappa_x(alpha, r, w0)*gamma_x(alpha, v**2, r, w0))


def gamma_x(s, v, r, w0, n=N):
    r = abs(r)

    if s == n-r+1:
        return gamma(1, v, n, r)

    return 1 + q_x(s, r, w0)*(v**2)/gamma_x(s+1, v, r, w0)


def gamma(l, v, n, r, n1=N1):
    """
        l belong {1, ..., n1}
        
    """
    r = abs(r)

    if l == n1 + 1:
        # v = x + iy
        x, y = v.real, v.imag
        # (1 + v**2)**0.5 = exp(i*(fi + fi_1)/2) * (1 + x**2 - y**2 + 4(x**2)(y**2))**0.5
        sqrt_expression = math.exp**((fi(v)+fi_1(v))*1j/2)*((1+x**2-y**2)**2 + 4*x**2*y**2)**0.5
        A = (-1 + sqrt_expression)/2
        return A + 1

    return 1 + a(l, v, n, r)/gamma(l+1, v, n, r)


def psi_x(j, v, r, w0):
    """
        j belong N
        abs(r) belong N0
    """
    result = 1

    for l in range(1, j+1):
        result *= h_x(l, v, r, w0)

    return result


def a(l, v, n, r):
    """
        l belongs {1, ..., n1}, n1 >> n
    """
    r = abs(r)
    return ((n+l)**2-r**2)*v**2/((2*(n+l)+1)*(2*(n+l+1)+1))
