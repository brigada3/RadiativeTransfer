#!/usr/bin/env python3
# -*-coding: utf-8 -*-

"""
     N   -- { 1, 2, 3 ... }
    N0  -- { 0, 1, 2, 3 ... }
"""

from sympy import *


N = 2001
N1 = 240*N
G = 0.999


s, r, w0, i, alpha, l, j = symbols('s r w0 i alpha l j')
f = G**i
epsilin_x = sqrt(s*(s+2*r))
kappa_x = (2*(s+abs(r))+1)*(1-w0*f.subs(i, s+r))


v = Symbol('v')
fi = Piecewise(
    (atan(im(v-1j)/re(v-1j)), re(v-1j) > 0),
    (-pi/2, And(Eq(re(v-1j), 0), im(v-1j) < 0)),
    (-pi/2 - atan(re(v-1j)/im(v-1j)), And(re(v-1j)<0, im(v-1j)<0)),
    (-pi, And(re(v-1j)<0, Eq(im(v-1j), 0))),
    (-pi + atan(im(v-1j)/re(v-1j)), And(re(v-1j)<0, im(v-1j)>0))
)

fi_1 = Piecewise(
    (atan(im(v+1j)/re(v+1j)), re(v+1j) > 0),
    (pi/2, And(Eq(re(v+1j), 0), im(v+1j) < 0)),
    (pi/2 - atan(re(v+1j)/im(v+1j)), And(re(v+1j)<0, im(v+1j)<0)),
    (pi, And(re(v+1j)<0, Eq(im(v+1j), 0))),
    (pi + atan(im(v+1j)/re(v+1j)), And(re(v+1j)<0, im(v+1j)>0))
)


q_x = Piecewise(
    ((s*(s+2*r))/((2*(s+r)+1)*(2*(s+r+1)+1)*(1-w0*f.subs(i, s+r)*(1-w0*f.subs(i, s+r+1)))), And(s+r+1<=N, N>=2)),
    ((s*(s+2*r))/((2*(s+r)+1)*(2*(s+r+1)+1)*(1-w0*f.subs(i, s+r))), Eq(s+r, N)),
    ((s*(s+2*r))/((2*(s+r)+1)*(2*(s+r+1)+1)), s >= N-r+1)
)

h_x = 1j*v*epsilon_x({s: alpha, r: r})/(kappa_x({s: alpha, r: r, w0: w0})*gamma_x({s: alpha, v: v**2, r: r, w0: w0}))
a = ((n+l)**2-r**2)*v**2/((2*(n+l)+1)*(2*(n+l+1)+1))

#TODO: Recursion!!!
gamma_x = Piecewise(
    (gamma.subs({l: 1, v: v, r: r}), Eq(s, n-r+1)),
    (1 + q_x.subs({s:s, r:r, w0:w0})*v/gamma_x.subs({s:s+1, v:v, r:r, w0:w0}), True)
)

gamma = Piecewise(
    ((-1+(exp(1)**(fi(v)+fi_1(v))*1j/2)*((1+re(v)**2-im(v)**2)**2 + 4*re(v)**2*im(v)**2)**0.5)/2 + 1, Eq(l, N1 + 1)),
    (1 + a.subs({l:l, v:v, r:r, w0:w0})/gamma.subs({l:l+1, v:v, r:r}), True)
)

psi_x = Product(h_x.subs({l:l, v:v, r:r, w0:w0}), (l, 1, j+1))


def P(j,r,m1):
    """
        j,r belong N0
    """
    r = abs(r)
    return ((factorial(j)/factorial(j)(j+2*r))**0.5)*P_legandr(j,r,m1)



#a belongs to {1,...,n-r},  n-r>=1
#var a,v,r,w0,m1

# def nu_x(a,v,r,w0,m1):
#     """
#     a belongs to {1,...,n-r},  n-r>=1
#     """
#     r = abs(r)
#     i = 1j
#     return (1/(kappa_x(a,r,w0)*gamma_x(a,v**2,r,w0)))*(i*v*epsilon_x(a+1,r)*nu_x(a+1,v,r,w0) + (2*(a+r) + 1)*fi(a+r)*P(a+r,r,m1))


# nu_x =  (1/(kappa_x(a,r,w0)*gamma_x(a,v**2,r,w0)))*(i*v*epsilon_x(a+1,r)*nu_x(a+1,v,r,w0) + (2*(a+r) + 1)*fi(a+r)*P(a+r,r,m1))

# gamma_x.subs({s:a, v:v**2, r:r, w0:w0})
# kappa_x.subs({s: a, r: r, w0: w0})
# epsilon_x.subs({s: a+1, r: r})
# nu_x.subs({a:a+1,v:v,r:r,w0:w0})
# fi.subs({v:a+r})
# P.subs({a:a+r,r:r,m1:m1}

#обратная рекурсия :(
#все символами вроде подставляется ( и i=1j, r=abs(r) )
nu_x =  (1/(kappa_x.subs({s: a, r: r, w0: w0})*gamma_x.subs({s:a, v:v**2, r:r, w0:w0})))*(i*v*epsilon_x.subs({s: a+1, r: r})*nu_x.subs({a:a+1,v:v,r:r,w0:w0}) + (2*(a+r) + 1)*fi.subs({v:a+r})*P.subs({a:a+r,r:r,m1:m1}))
