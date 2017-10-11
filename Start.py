#!/usr/bin/env python3
# -*-coding: utf-8 -*-
import math

f = []
n = len(f)

def epsilon(s, r):
    r = abs(r)
    return (s*(s+2*r))**0.5

def kappa(s, r, w):
    r = abs(r)
    return(2*(s+r)+1)*(1-w*f[s+r])

def fi(nu):
    nu = nu - 1j
    if nu.real > 0:
        return math.atan(nu.imag/nu.real)
    elif nu.real == 0 and nu.imag < 0:
        return -math.pi/2
    elif nu.real < 0 and nu.imag < 0:
        return -math.pi/2 - math.atan(nu.real/nu.imag)
    elif nu.real < 0 and nu.imag == 0:
        return -math.pi
    elif nu.real < 0 and nu.imag > 0:
        return -math.pi + math.atan(nu.imag/nu.real)


def fi_1(nu):
    nu = nu + 1j
    if nu.real > 0:
        return math.atan(nu.imag/nu.real)
    elif nu.real == 0 and nu.imag < 0:
        return math.pi/2
    elif nu.real < 0 and nu.imag < 0:
        return math.pi/2 - math.atan(nu.real/nu.imag)
    elif nu.real < 0 and nu.imag == 0:
        return math.pi
    elif nu.real < 0 and nu.imag > 0:
        return math.pi + math.atan(nu.imag/nu.real)


def q(s,r,w):
    r = abs(r)
    result = s*(s+2*r)
    if s+1+r <= n and n >= 2:
        return result/((2*(s+r)+1)*(2*(s+r+1)+1)*(1-w*f[s+r])*(1-w*f[s+r+1]))
    elif s+r == n:
        return result/((2*(s+r)+1)*(2*(s+r+1)+1)*(1-w*f[s+r]))
    elif s >= n-r+1:
        return result/((2*(s+r)+1)*(2*(s+r+1)+1))

def h(a,v,r,w):
    r = abs(r)
    return 1j*v*epsilon(a,r)/(kappa(a,r,w)*gamma(a,v**2,r,w))

def gamma(a,v,r,w):
    pass

def psi(j,v,r,w):
    result = 1;
    for l in range(1, j+1):
        result *= h(l, v, r, w)

    return result
