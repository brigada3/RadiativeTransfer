#!/usr/bin/env python3
# -*-coding: utf-8 -*-

"""
     N   -- { 1, 2, 3 ... }
    N0  -- { 0, 1, 2, 3 ... }
"""

import sympy 


N = 2001
N1 = 240*N
G = 0.999


s, r, w0, i, n = sympy.symbols('s r w0 i n')
alpha, l, j, v = sympy.symbols('alpha l j v')

f = G**i
epsilon_x = sympy.sqrt(s*(s+2*r))
kappa_x = (2*(s+r)+1)*(1-w0*f.subs(i, s+r))

fi = sympy.Piecewise(
    (sympy.atan(sympy.im(v-1j)/sympy.re(v-1j)), sympy.re(v-1j) > 0),
    (-sympy.pi/2, sympy.And(sympy.Eq(sympy.re(v-1j), 0), sympy.im(v-1j) < 0)),
    (-sympy.pi/2 - sympy.atan(sympy.re(v-1j)/sympy.im(v-1j)), sympy.And(sympy.re(v-1j)<0, sympy.im(v-1j)<0)),
    (-sympy.pi, sympy.And(sympy.re(v-1j)<0, sympy.Eq(sympy.im(v-1j), 0))),
    (-sympy.pi + sympy.atan(sympy.im(v-1j)/sympy.re(v-1j)), sympy.And(sympy.re(v-1j)<0, sympy.im(v-1j)>0))
)

fi_1 = sympy.Piecewise(
    (sympy.atan(sympy.im(v+1j)/sympy.re(v+1j)), sympy.re(v+1j) > 0),
    (sympy.pi/2, sympy.And(sympy.Eq(sympy.re(v+1j), 0), sympy.im(v+1j) < 0)),
    (sympy.pi/2 - sympy.atan(sympy.re(v+1j)/sympy.im(v+1j)), sympy.And(sympy.re(v+1j)<0, sympy.im(v+1j)<0)),
    (sympy.pi, sympy.And(sympy.re(v+1j)<0, sympy.Eq(sympy.im(v+1j), 0))),
    (sympy.pi + sympy.atan(sympy.im(v+1j)/sympy.re(v+1j)), sympy.And(sympy.re(v+1j)<0, sympy.im(v+1j)>0))
)

q_x = sympy.Piecewise(
    ((s*(s+2*r))/((2*(s+r)+1)*(2*(s+r+1)+1)*(1-w0*f.subs(i, s+r)*(1-w0*f.subs(i, s+r+1)))), sympy.And(s+r+1<=N, N>=2)),
    ((s*(s+2*r))/((2*(s+r)+1)*(2*(s+r+1)+1)*(1-w0*f.subs(i, s+r))), sympy.Eq(s+r, N)),
    ((s*(s+2*r))/((2*(s+r)+1)*(2*(s+r+1)+1)), s >= N-r+1)
)

def get_gamma(l=l, v=v, n=n, r=r):
    if l.is_Number and v.is_Number and n.is_Number and r.is_Number:
        if l == N1 + 1:
            return (-1+(exp(1)**(fi.subs({v: v})+fi_1.subs({v: v}))*1j/2)*((1+sympy.re(v)**2-sympy.im(v)**2)**2 + 4*sympy.re(v)**2*sympy.im(v)**2)**0.5)/2 + 1
        else:
            return 1+a.subs( {l:l, v:v, r:r, w0:w0, n:N} )/get_gamma(l+1, v, n, r).subs( {l:l+1, v:v, r:r, n:n})
    else:
        return 1+a.subs( {l:l, v:v, r:r, w0:w0, n:N} )/sympy.Function('gamma')(l+1, v, n, r)
        
def get_gamma_x(s=s, v=v, r=r, w0=w0):
    if s.is_Number and v.is_Number and r.is_Number and w0.is_Number:
        if s == N-r+1:
            return gamma.subs({l: 1, v: v, r: r})
        else:
            return 1 + q_x.subs({s:s, r:r, w0:w0})*v/get_gamma_x(s+1, v, r, w0).subs({s: s+1, v: v, r: r, w0: w0})
    else:
        return 1 + q_x.subs({s:s, r:r, w0:w0})*v/sympy.Function('gamma_x')(s+1, v, r, w0)

h_x = (1j*v*epsilon_x.subs( {s: alpha, r: r } )/
    (kappa_x.subs( {s: alpha, r: r, w0: w0} )*
    get_gamma_x(alpha, v**2, r, w0)).subs( {s: alpha, v: v**2, r: r, w0: w0}))

a = ((n+l)**2-r**2)*v**2/((2*(n+l)+1)*(2*(n+l+1)+1))

psi_x = sympy.Product(h_x.subs({l:l, v:v, r:r, w0:w0}), (l, 1, j+1))
