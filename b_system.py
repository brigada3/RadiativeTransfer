#!/usr/bin/env python3
# -*-coding: utf-8 -*-

import sympy


N = 2001
N1 = 240*N
G = 0.999


s, r, w0, i, n = sympy.symbols('s r w0 i n')
alpha, l, j, v = sympy.symbols('alpha l j v')
mu = sympy.symbols('mu')


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
    if l == N1 + 1:
        return (
        (-1 + (exp(1)**(fi.subs({v: v})+fi_1.subs({v: v}))*1j/2) *
        ( (1+sympy.re(v)**2-sympy.im(v)**2)**2 + 4*sympy.re(v)**2*sympy.im(v)**2)**0.5  ) / 2 + 1
    )
    elif l.is_Number:
        gamma_l_plus_1 = get_gamma(l+1, v, n, r).subs( {l:l+1, v:v, r:r, n:n})
    else:
        gamma_l_plus_1 = sympy.Function('gamma')(l+1, v, n, r)

    return 1+a.subs( {l:l, v:v, r:r, w0:w0, n:N} ) / gamma_l_plus_1


def get_gamma_x(s=s, v=v, r=r, w0=w0):
    if s == N-r+1:
        return gamma.subs( {l: 1, v: v, r: r} )
    elif s.is_Number and r.is_Number:
        gamma_x_s_plus_1 = get_gamma_x(s+1, v, r, w0).subs( {s: s+1, v: v, r: r, w0: w0} )
    else:
        gamma_x_s_plus_1 = sympy.Function('gamma_x')(s+1, v, r, w0)

    return 1 + q_x.subs( {s:s, r:r, w0:w0} )*v/gamma_x_s_plus_1


h_x = (1j*v*epsilon_x.subs( {s: alpha, r: r } )/
    (kappa_x.subs( {s: alpha, r: r, w0: w0} )*
    get_gamma_x(alpha, v**2, r, w0)).subs( {s: alpha, v: v**2, r: r, w0: w0}))

a = ((n+l)**2-r**2)*v**2/((2*(n+l)+1)*(2*(n+l+1)+1))


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



psi_x = sympy.Product(h_x.subs({l:l, v:v, r:r, w0:w0}), (l, 1, j+1))

def get_Px(i=s+r+1, r=r, mu=mu):
    if i == r-1:
        return 0
    elif i == r:
        return (1-mu**2)**(r/2) * (2)**(-r) * sympy.sqrt(sympy.factorial(2*r)/sympy.factorial(r)**2)
    elif i.is_Number and r.is_Number:
        Px_i_minus_1 = get_Px(i-1, r, mu).subs( {i: i-1, r: r, mu: mu} )
        Px_i_minus_2 = get_Px(i-2, r, mu).subs( {i: i-2, r: r, mu: mu} )
    else:
        Px_i_minus_1 = sympy.Function('Px')(i-1, r, mu)
        Px_i_minus_2 = sympy.Function('Px')(i-2, r, mu)

    return (
        ( 2*i*mu*Px_i_minus_1 - sympy.sqrt((i-r-1)*(i+r-1)) * Px_i_minus_2 ) *
        ( (i-r)*(i+r) )**(-1/2)
    )

